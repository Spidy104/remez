$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path $PSScriptRoot -Parent
$originalPath = $env:PATH
$originalUbsan = $env:UBSAN_OPTIONS
$originalAsan = $env:ASAN_OPTIONS

function Invoke-Xmake {
    & xmake @args
    if ($LASTEXITCODE -ne 0) { throw "xmake failed with exit $LASTEXITCODE" }
}

function Test-Probe($executable, $mode, $expected) {
    $report = (& $executable $mode 2>&1 | Out-String)
    $exitCode = $LASTEXITCODE
    $report | Set-Content -LiteralPath (Join-Path (Split-Path $executable) "$mode.log")
    if ($mode -eq 'clean') {
        if ($exitCode -ne 0 -or $report.Trim() -ne '42') { throw "Clean probe failed: $report" }
    } elseif ($exitCode -eq 0 -or $report -notmatch $expected) {
        throw "Expected $expected for $mode; exit=$exitCode`n$report"
    }
    Write-Host "PASS $mode"
}

Push-Location -LiteralPath $projectRoot
try {
    $vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
    $vsRoot = (& $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath).Trim()
    $toolVersion = (Get-Content -LiteralPath (Join-Path $vsRoot 'VC\Auxiliary\Build\Microsoft.VCToolsVersion.default.txt')).Trim()
    $msvcRuntime = Join-Path $vsRoot "VC\Tools\MSVC\$toolVersion\bin\HostX64\x64"
    $clangResource = (& clang++ --print-resource-dir).Trim()
    if ($LASTEXITCODE -ne 0) { throw 'Cannot locate Clang runtime' }
    $env:UBSAN_OPTIONS = 'halt_on_error=1:print_stacktrace=1'
    $env:ASAN_OPTIONS = 'halt_on_error=1'

    foreach ($compiler in @('msvc', 'clang')) {
        $policy = if ($compiler -eq 'clang') { 'build.sanitizer.address,build.sanitizer.undefined' } else { 'build.sanitizer.address' }
        $runtime = if ($compiler -eq 'clang') { Join-Path $clangResource 'lib\windows' } else { $msvcRuntime }
        $env:PATH = "$runtime;$originalPath"
        Invoke-Xmake f -c -y -m debug "--toolchain=$compiler" "--builddir=build/$compiler" "--policies=$policy"
        Invoke-Xmake -y sanitizer_probe
        Invoke-Xmake -y tests
        Invoke-Xmake run tests
        $executable = Join-Path $projectRoot "build\$compiler\windows\x64\debug\sanitizer_probe.exe"
        Write-Host "$compiler sanitizer verification"
        Test-Probe $executable 'clean' ''
        Test-Probe $executable 'bounds' 'AddressSanitizer: heap-buffer-overflow'
        Test-Probe $executable 'use-after-free' 'AddressSanitizer: heap-use-after-free'
        if ($compiler -eq 'clang') {
            Test-Probe $executable 'overflow' 'runtime error: signed integer overflow'
        }
    }
    $env:PATH = "$msvcRuntime;$originalPath"
    Invoke-Xmake f -c -y -m debug --toolchain=msvc --builddir=build/msvc-crt --memory_checks=y
    Invoke-Xmake -y memory_tests
    Invoke-Xmake run memory_tests
    foreach ($compiler in @('msvc', 'clang')) {
        foreach ($mode in @('debug', 'release')) {
            Invoke-Xmake f -c -y "-m$mode" "--toolchain=$compiler" "--builddir=build/$compiler-$mode"
            Invoke-Xmake -y tests
            Invoke-Xmake run tests
        }
    }
} finally {
    $env:PATH = $originalPath
    $env:UBSAN_OPTIONS = $originalUbsan
    $env:ASAN_OPTIONS = $originalAsan
    Pop-Location
}




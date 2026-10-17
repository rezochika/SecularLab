<#
.SYNOPSIS
Build project-local GMP C++/Boost dependencies using an existing CLion MinGW.
.DESCRIPTION
Requires Git for Windows. Downloads verified upstream source archives to
build/deps, builds with the existing compiler, and creates local CMake presets
when CMakeUserPresets.json does not already exist. Installs no system packages.
#>
[CmdletBinding()]
param(
    [string]$ClionRoot,
    [string]$GitBash,
    [switch]$NoPresets
)

$ErrorActionPreference = 'Stop'
if ([Environment]::OSVersion.Platform -ne [PlatformID]::Win32NT) {
    throw 'Run this entrypoint in Windows PowerShell or PowerShell on Windows.'
}
$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path

function Test-ClionRoot([string]$candidate) {
    foreach ($relative in @('bin/mingw/bin/g++.exe', 'bin/mingw/include/gmp.h',
            'bin/mingw/lib/libgmp.a', 'bin/cmake/win/x64/bin/cmake.exe',
            'bin/ninja/win/x64/ninja.exe')) {
        if (-not (Test-Path -LiteralPath (Join-Path $candidate $relative) -PathType Leaf)) { return $false }
    }
    return $true
}

if (-not $ClionRoot) {
    $jetbrainsRoot = Join-Path $env:ProgramFiles 'JetBrains'
    $candidates = @()
    if (Test-Path -LiteralPath $jetbrainsRoot -PathType Container) {
        $candidates = @(Get-ChildItem -LiteralPath $jetbrainsRoot -Directory -Filter 'CLion*' |
            Sort-Object LastWriteTime -Descending)
    }
    foreach ($candidate in $candidates) {
        if (Test-ClionRoot $candidate.FullName) { $ClionRoot = $candidate.FullName; break }
    }
    if (-not $ClionRoot) {
        throw 'CLion with bundled MinGW/GMP was not found. Pass -ClionRoot with its installation folder.'
    }
}
$ClionRoot = (Resolve-Path -LiteralPath $ClionRoot).Path
if (-not (Test-ClionRoot $ClionRoot)) {
    throw 'The selected CLion folder must contain bundled MinGW GCC, static GMP, CMake, and Ninja.'
}

if (-not $GitBash) {
    $gitCandidates = @((Join-Path $env:ProgramFiles 'Git/bin/bash.exe'))
    $gitCommand = Get-Command git.exe -ErrorAction SilentlyContinue
    if ($gitCommand) {
        $gitCandidates += Join-Path (Split-Path (Split-Path $gitCommand.Source -Parent) -Parent) 'bin/bash.exe'
    }
    foreach ($candidate in $gitCandidates) {
        if (Test-Path -LiteralPath $candidate -PathType Leaf) { $GitBash = $candidate; break }
    }
    if (-not $GitBash) { throw 'Git for Windows Bash was not found. Pass -GitBash with the existing bash.exe path.' }
}
$GitBash = (Resolve-Path -LiteralPath $GitBash).Path
$shellScript = Join-Path $PSScriptRoot 'setup-clion-deps.sh'
Write-Host "Using CLion: $ClionRoot"
& $GitBash $shellScript $ClionRoot
if ($LASTEXITCODE -ne 0) { throw "Dependency setup failed (exit $LASTEXITCODE). See build/deps for logs." }

$presetsPath = Join-Path $projectRoot 'CMakeUserPresets.json'
if (-not $NoPresets -and -not (Test-Path -LiteralPath $presetsPath)) {
    $clion = $ClionRoot.Replace('\', '/')
    $release = [ordered]@{
        name = 'windows-clion'
        displayName = 'Windows - CLion bundled MinGW Release'
        generator = 'Ninja'
        binaryDir = '${sourceDir}/build/windows-clion'
        condition = @{ type = 'equals'; lhs = '${hostSystemName}'; rhs = 'Windows' }
        cacheVariables = [ordered]@{
            CMAKE_BUILD_TYPE = 'Release'
            BUILD_TESTING = 'ON'
            CMAKE_CXX_COMPILER = "$clion/bin/mingw/bin/g++.exe"
            CMAKE_MAKE_PROGRAM = "$clion/bin/ninja/win/x64/ninja.exe"
            CMAKE_PREFIX_PATH = '${sourceDir}/build/deps/mingw'
        }
        environment = @{ PATH = "$clion/bin/mingw/bin;" + '$penv{PATH}' }
        vendor = @{ 'jetbrains.com/clion' = @{ toolchain = 'MinGW' } }
    }
    $debugPreset = [ordered]@{
        name = 'windows-clion-debug'
        displayName = 'Windows - CLion bundled MinGW Debug'
        inherits = 'windows-clion'
        binaryDir = '${sourceDir}/build/windows-clion-debug'
        cacheVariables = @{ CMAKE_BUILD_TYPE = 'Debug' }
    }
    $presets = [ordered]@{
        version = 3
        configurePresets = @($release, $debugPreset)
        buildPresets = @(
            @{ name = 'windows-clion'; configurePreset = 'windows-clion' },
            @{ name = 'windows-clion-debug'; configurePreset = 'windows-clion-debug' }
        )
        testPresets = @(
            @{ name = 'windows-clion'; configurePreset = 'windows-clion'; output = @{ outputOnFailure = $true } },
            @{ name = 'windows-clion-debug'; configurePreset = 'windows-clion-debug'; output = @{ outputOnFailure = $true } }
        )
    }
    $json = $presets | ConvertTo-Json -Depth 10
    [IO.File]::WriteAllText($presetsPath, $json + [Environment]::NewLine, [Text.UTF8Encoding]::new($false))
    Write-Host 'Created local windows-clion and windows-clion-debug presets.'
} elseif (-not $NoPresets) {
    Write-Host 'Existing CMakeUserPresets.json retained.'
}
Write-Host 'Dependencies are ready. See docs/WINDOWS_SETUP.md for CLion and command-line build steps.'

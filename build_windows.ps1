$ErrorActionPreference = "Stop"

$vsVarsPath = "${env:ProgramFiles}\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"

if (-Not (Test-Path $vsVarsPath)) {
    Write-Error "vcvars64.bat not found. Make sure you have Visual Studio with C++ tools installed."
    exit 1
}

cmd /c "`"$vsVarsPath`" && cl /W4 /DFOURIER_IMPLEMENTATION main.c"

if ($args -contains "run") {
    Write-Output "The 'run' argument was provided. Executing command..."
    cmd /c "main.exe"
}

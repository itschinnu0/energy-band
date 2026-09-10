# Energy Band Diagram Simulator

## Deployment and Packaging

This project uses PyInstaller to build a Windows executable.

### Packaging Command

To build the executable, run:

`powershell
uv run pyinstaller --noconfirm --onedir --clean --copy-metadata streamlit --add-data "app.py;." launcher.py
`

The resulting executable will be available at dist/launcher/launcher.exe.

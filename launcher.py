"""Packaging/runtime launcher."""
import os
import sys
import streamlit.web.cli as stcli

def resolve_path(path):
    if getattr(sys, 'frozen', False):
        return os.path.join(sys._MEIPASS, path)
    return os.path.join(os.path.dirname(__file__), path)

def main():
    app_path = resolve_path("app.py")
    sys.argv = [
        "streamlit",
        "run",
        app_path,
        "--global.developmentMode=false"
    ]
    sys.exit(stcli.main())

if __name__ == "__main__":
    main()

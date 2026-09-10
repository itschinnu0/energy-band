"""Packaging/runtime launcher."""
import os
import sys

# Explicit imports to ensure PyInstaller bundles all Streamlit runtime components
import streamlit
import streamlit.runtime.scriptrunner.magic_funcs
import streamlit.runtime.scriptrunner.magic
import streamlit.runtime.scriptrunner.script_runner
import streamlit.runtime.scriptrunner.exec_code
import streamlit.runtime.scriptrunner_utils.script_run_context
import streamlit.web.cli as stcli


def resolve_path(path):
    if getattr(sys, "frozen", False):
        return os.path.join(sys._MEIPASS, path)
    return os.path.join(os.path.dirname(__file__), path)


def main():
    app_path = resolve_path("app.py")
    sys.argv = [
        "streamlit",
        "run",
        app_path,
        "--global.developmentMode=false",
    ]
    sys.exit(stcli.main())


if __name__ == "__main__":
    main()

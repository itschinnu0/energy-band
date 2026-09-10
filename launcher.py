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


# Explicit imports to ensure PyInstaller bundles all local project packages and Streamlit runtime components
import physics
import visualization
import app


def main():
    # When frozen, ensure the extraction root is at the front of sys.path so subprocesses/exec discover local packages
    base_dir = resolve_path(".")
    if base_dir not in sys.path:
        sys.path.insert(0, base_dir)

    app_path = resolve_path("app.py")
    # Pass through user-provided command-line arguments (e.g., --server.port, --server.headless)
    extra_args = sys.argv[1:]
    sys.argv = [
        "streamlit",
        "run",
        app_path,
        "--global.developmentMode=false",
    ] + extra_args
    sys.exit(stcli.main())


if __name__ == "__main__":
    main()

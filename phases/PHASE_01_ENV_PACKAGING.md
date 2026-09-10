# Phase 1 — Environment and Packaging Feasibility Spike

## Objective
Prove that the chosen Python + Streamlit + PyInstaller deployment path works before substantial physics code is written.

## Tasks
1. Initialize/manage dependencies with `uv`.
2. Pin Python 3.12.
3. Add:
   - numpy
   - scipy
   - matplotlib
   - streamlit
   - pytest
4. Create a minimal Streamlit page.
5. Create `launcher.py` that starts the local Streamlit server and opens the browser.
6. Build a Windows executable using PyInstaller.
7. Record hidden imports/assets required for the frozen application.
8. Keep this spike minimal.

## Important
Do not solve physics in this phase.

## Acceptance
- Development server launches.
- Minimal page renders.
- Frozen executable starts the application.
- Browser auto-launch works.
- Packaging instructions are documented.
- Any PyInstaller issue is documented rather than patched with unexplained hacks.

## Suggested commands
- `uv sync`
- `uv run pytest`
- `uv run streamlit run app.py`
- Windows packaging command documented in the project README.

## Deliverable
A minimal but independently launchable application and reproducible packaging instructions.

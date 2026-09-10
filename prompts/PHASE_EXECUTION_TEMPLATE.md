# Phase Execution Prompt

Use this prompt with a coding agent after selecting a phase.

> Read `AGENTS.md`, `IMPLEMENTATION_PLAN_SOURCE.md`, and the selected phase file under `phases/`.
>
> You are implementing this phase only.
>
> First inspect the repository and identify existing work. Then implement the phase with production-quality Python, explicit units, type hints, tests, and clear error handling.
>
> Do not invent unsupported physics. Do not silently clip invalid states. Do not weaken tests to force a pass.
>
> At the end:
> 1. run the phase tests;
> 2. run the full suite if practical;
> 3. report files changed;
> 4. report physics equations implemented;
> 5. report tests and results;
> 6. report unresolved issues;
> 7. state whether the phase gate is passed.
>
> If the specification is physically ambiguous or inconsistent, stop the affected implementation and explain the exact issue rather than guessing.

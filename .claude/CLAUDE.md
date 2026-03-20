# Claude Instructions for E3SM Diagnostics

> **Canonical rules:** [`AGENTS.md`](../AGENTS.md) at the repository root.
> All coding standards, architecture constraints, testing philosophy, dependency
> policies, and concurrency rules are defined there. Follow them in all
> responses.

## Claude-Specific Guidance

- Always apply the rules in `AGENTS.md` when generating or reviewing code.
- Prefer patterns already present in the repository over speculative or
  generic approaches.
- Do not generate references to files, modules, configuration keys, or
  dependencies that do not exist in the repository.
- When uncertain about a convention, consult `AGENTS.md` or existing source
  code rather than guessing.
- Do not introduce new dependencies, tools, or architectural patterns unless
  explicitly requested and justified.

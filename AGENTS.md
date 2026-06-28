# alphafold3


<!-- AGENTS-CODEX-PREFLIGHT-v1 -->
## For Codex / non-Claude agents — read this repo's CLAUDE.md, then prove context

If you are Codex (or any agent that does not auto-load `CLAUDE.md`):
- **Read `./CLAUDE.md` in THIS directory, in full, before any work.** It carries this repo's
  binding invariants, the contracts that apply, gotchas, and build/test commands. This
  `AGENTS.md` is only a pointer (summary + links) — the repo-specific rules live in `./CLAUDE.md`
  and in `~/repos/governance/`.
- **Run the mandatory PHASE 0.25 context proof when you start a task in a new repo:**
  `python3 ~/repos/governance/enforcement/generate_challenge.py`, emit a `CONTEXT_PROOF`
  (format: `~/repos/governance/enforcement/test_fixtures/complete_proof.yaml`), then
  `python3 ~/repos/governance/enforcement/verify_agent_context.py --proof <file>` — proceed
  only on exit 0. Authoritative: `governance/protocols/MANDATORY_CONTEXT_VERIFICATION_PROTOCOL.md`.
- Treat `~/repos/governance/INDEX.md` and its contracts/policies as binding.
<!-- /AGENTS-CODEX-PREFLIGHT-v1 -->


<!-- GOVERNANCE-PREFLIGHT-v1 -->
## Governance Pre-Flight (summary — binding rules live in governance/)

All agents — Claude, Codex, Grok, Gemini, Hermes — before starting a task:
- Complete the startup audit and the **PHASE 0.5 pre-flight restatement**: restate to the
  user a 3–5 step plan plus the three most relevant governance policies, before doing the work.
- Use the **canonical document template** for any document — do not invent a format.

This is a summary; the binding rules and full checklists live in governance (source of truth):
- `~/repos/governance/policies/AGENT_INTERACTION_POLICY.md` — startup sequence + PHASE 0.5
- `~/repos/governance/standards/DOCUMENT_TEMPLATE_REGISTRY.md` — which template to use
- `~/repos/governance/INDEX.md` — master registry of all contracts, policies, gates
<!-- /GOVERNANCE-PREFLIGHT-v1 -->

Guidance for agents working in this repository.


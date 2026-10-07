# Operations: running the agents long-term

This guide is for the owner. The agents follow CLAUDE.md and [agent-charter.md](agent-charter.md), and this page
shows how to keep them running without re-pasting the charter. Everything here uses built-in Claude Code features
(project skills, `/loop`, `/schedule` routines, Desktop scheduled tasks, `/batch`, `claude-code-action`) rather than
custom orchestration code. The reasons are in DECISIONS.md D-009. All behaviour described here was checked against
the official docs on 2026-10-07 (links in SOURCES.md, "Agent tooling").

| Skill | What it does |
|---|---|
| `/research-cycle [focus] [--dry-run]` | One complete work cycle: orient, review replies, WIP cap, choose action and mode, implement, test, PR, handoff |
| `/reuse-check <need>` | Reuse-vs-build search in an isolated subagent. Returns a D-NNN entry and SOURCES lines |
| `/vet-candidate <source_uid> [run_id]` | Tests ordinary explanations for a candidate and writes a vetting record |

The skill files live in `.claude/skills/`. Because they're committed, they load locally, in cloud routines and in
GitHub Actions. `.claude/settings.json` holds the project permission rules (section 8).

**Trust the repository once.** Start `claude.exe` interactively in the repository root and accept the workspace trust
dialog. The `allow` rules in `.claude/settings.json` apply only after that, and the trust also covers worktrees of
the checkout. `-p` and SDK runs never show the dialog, and that includes the GitHub Action, so they ignore project
`allow` rules in an untrusted folder. That's why `research-cycle` lists the cycle's tools in its own `allowed-tools`
(skill grants are not trust-gated) and `claude.yml` passes `--allowedTools`. `deny` rules apply everywhere.

**Only your words are instructions.** The skills treat issues, comments and PR text by anyone other than you as
data. So steer the agents through TASKS.md, your own issues and your own review comments.

## 1. One cycle

```text
/research-cycle                      # the agent picks the highest-value action
/research-cycle M1 integration       # with a focus
/research-cycle --dry-run            # orientation plus the choice it would make; changes nothing
```

Headless, from Git Bash: `claude.exe -p "/research-cycle --dry-run" --max-turns 15 --permission-mode plan`. The
2026-10-07 test run finished in 5 turns.

## 2. Local long run: `/loop`

`/loop` re-runs a prompt inside an open session.

```text
/loop /research-cycle          # self-paced: Claude picks 1 min to 1 h between iterations
/loop 2h /research-cycle       # fixed interval
```

Limits (see [scheduled tasks](https://code.claude.com/docs/en/scheduled-tasks)):
- It runs only while the session is open and idle. Missed fires are not caught up.
- Recurring tasks expire after 7 days.
- `--resume` restores fixed-interval tasks but not a self-paced loop.
- `Esc` stops a self-paced loop. Ask Claude to cancel a fixed one.
- All iterations share one context. Each cycle externalizes its knowledge (skill step 9), but for runs longer than a
  day prefer a cloud routine, which starts a fresh session every run.
- A scheduled fire runs only skills that Claude may invoke itself. That is why `research-cycle` must not set
  `disable-model-invocation`.
- Local runs stop at permission prompts for commands outside the allow list (section 8). Approve them with
  "don't ask again", or add the rule to `.claude/settings.local.json`.

## 3. Cloud routine via `/schedule` (recommended primary scheduler)

[Routines](https://code.claude.com/docs/en/routines) run on Anthropic's cloud with your laptop off. Each run is a
fresh clone of the repository in a fresh session. Routines are a research preview, so details may change.

1. **Connect GitHub.** Run `/web-setup` once, or install the [Claude GitHub App](https://github.com/apps/claude).
2. **Create a cloud environment** at [claude.ai/code](https://claude.ai/code) (environment selector, then
   **Add cloud environment**). Configure it like this:
   - **Network access: Custom.** Check **Also include default list of common package managers**. The default
     list already covers PyPI, GitHub and `*.amazonaws.com`, which includes the STScI S3 bucket
     `stpubdata.s3.amazonaws.com`. Then add these allowed domains:
     ```text
     *.stsci.edu
     *.cds.unistra.fr
     cdsxmatch.u-strasbg.fr
     ned.ipac.caltech.edu
     gea.esac.esa.int
     arxiv.org
     export.arxiv.org
     raw.githubusercontent.com
     s3.amazonaws.com
     doi.org
     zenodo.org
     jwst-uncover.github.io
     www.fe.infn.it
     ```
     These are the hosts we observed (2026-10-07, astroquery 0.4.11, by logging every DNS lookup) for a MAST query,
     a product list and a `_cat.ecsv` download (`mast.stsci.edu`), an S3 byte-range FITS read
     (`stpubdata.s3.amazonaws.com`), SIMBAD (`simbad.cds.unistra.fr`), VizieR (`vizier.cds.unistra.fr`), CDS XMatch
     (`cdsxmatch.u-strasbg.fr`), NED (`ned.ipac.caltech.edu`) and the Gaia archive (`gea.esac.esa.int`). The
     wildcards also cover the other MAST hosts in astroquery's `mast` module (`catalogs.mast.stsci.edu`,
     `auth.mast.stsci.edu`) and the alternative XMatch host `cdsxmatch.cds.unistra.fr` (also live). The rest serve
     literature checks (arXiv, DOI links), raw files from public repositories (the pinned Mahler+2022 model), the
     path-style DJA catalog URLs (`s3.amazonaws.com/grizli-v2/...`) and the published lens models of the cluster
     fields (SOURCES.md: Scofield+2025 on Zenodo, UNCOVER, Bergamini+2023). `raw.githubusercontent.com` and
     `s3.amazonaws.com` are listed explicitly: harmless if the default list already covers them, and needed if it
     does not. The UNCOVER maps sit behind a Google Drive link; add `drive.google.com` and
     `drive.usercontent.google.com` only if that file is needed. A blocked request fails with `403` and
     `x-deny-reason: host_not_allowed`. Add the host when a run reports one, for example a docs site such as
     `*.readthedocs.io`.
   - **Environment variables:** none. Never put secrets here, because anyone using the environment can read them.
     A MAST token isn't needed for public data.
   - **Setup script:** optional. uv, Python 3 and `gh` are preinstalled, and the cycle creates `.venv` itself. To warm
     uv's cache, use:
     ```bash
     #!/bin/bash
     # Cached about 7 days. It must exit 0 and finish in about 5 minutes.
     command -v uv >/dev/null 2>&1 || python3 -m pip install --quiet uv || true
     uv venv /opt/warm-venv --python python3 --quiet && \
       uv pip install --python /opt/warm-venv/bin/python --quiet \
         astropy astroquery numpy scipy pandas pyarrow scikit-learn matplotlib pyyaml pytest ruff fsspec s3fs || true
     exit 0
     ```
3. **Create the routine.** Run `/schedule` in a local session (or open [claude.ai/code/routines](https://claude.ai/code/routines)).
   Choose the repository `DD-Ching/jwst-anomaly-research`, the environment from step 2 and no connectors (new
   routines attach every connected connector; remove them). The plain prompt `/research-cycle` works; the live
   routine (since 2026-10-08: hourly at :07 UTC, Opus 5.5) uses the longer prompt in
   [cloud-routine-prompt.md](cloud-routine-prompt.md), which loops cycles for about 40 minutes and adds the
   D-023 mission and the cloud-specific rules. Start slower (daily, then cron `7 */8 * * *`) if your review can't keep up. The minimum
   interval is 1 hour; start a few minutes past the hour, since on-the-hour starts can lag. The WIP cap (3 open agent
   PRs) limits unmerged work. Runs end after about 40 minutes, so consecutive hourly runs don't overlap; a run
   leaves alone PRs labelled `local-wip` and branches with a commit in the last 15 minutes.
4. **Check the first run** with **Run now**, then open the session. A green status only means the session exited
   cleanly. Read the transcript, or ask `/schedule why did my research cycle do nothing?`.

**GitHub from a cloud session.** GraphQL is blocked there (`HTTP 403: GitHub GraphQL is not available from Claude
Code sessions`), so every `gh pr ...` and `gh issue ...` command fails. Two routes work:
- the session's GitHub MCP tools (load them with ToolSearch: `mcp__github__create_pull_request`,
  `pull_request_read`, `merge_pull_request`). The first routine run opened #38 with them, without a permission
  prompt;
- REST through `gh api`, with literal paths (shell variables don't persist between tool calls):

| Purpose | Cloud call |
|---|---|
| open PRs | `gh api 'repos/DD-Ching/jwst-anomaly-research/pulls?state=open'` |
| one PR's state | `gh api repos/DD-Ching/jwst-anomaly-research/pulls/N --jq '[.user.login, .head.repo.full_name, .head.ref, .head.sha, .mergeable_state] \| @tsv'` |
| a PR's files (every page) | `gh api --paginate repos/DD-Ching/jwst-anomaly-research/pulls/N/files --jq '.[].filename'` |
| open issues | `gh api 'repos/DD-Ching/jwst-anomaly-research/issues?state=open'` (entries with `pull_request` are PRs) |
| comments and reviews | `gh api repos/DD-Ching/jwst-anomaly-research/issues/N/comments`, `.../pulls/N/reviews`, `.../pulls/N/comments` |
| comment on a PR or issue | `gh api repos/DD-Ching/jwst-anomaly-research/issues/N/comments -f body="..."` |
| reply to an inline comment | `gh api repos/DD-Ching/jwst-anomaly-research/pulls/N/comments/ID/replies -f body="..."` |
| open an issue | `gh api repos/DD-Ching/jwst-anomaly-research/issues -f title="..." -f body="..."` |
| open a PR | `mcp__github__create_pull_request`, or `gh api repos/DD-Ching/jwst-anomaly-research/pulls -f title="..." -f head=claude/<slug> -f base=main -f body="..."` |
| add labels | `gh api repos/DD-Ching/jwst-anomaly-research/issues/N/labels -f "labels[]=agent"` (adds; MCP `issue_write` replaces the whole label set, so don't label with it) |
| wait for CI | the loop below, as one background Bash call |
| merge | `mcp__github__merge_pull_request` with `merge_method: squash` and `expectedHeadSha` = the reviewed head (the REST merge API is denied by `.claude/settings.json`) |

Write each body in the command itself (`-f body="$(cat <<'EOF'` on one line, the text, then `EOF` and `)"` on
their own lines); `=@file` and `--input` are denied so local files can't be posted by accident.

**Conflicts.** A PR that conflicts with `main` has `mergeable_state` `dirty` and gets no `pull_request` CI at all
(#38's first three commits had none). Merge `origin/main` into the branch, review the conflict resolution
(`git diff HEAD^1 HEAD` on the conflicted files), test and push. The merge is then pinned to the new head.

**Waiting for CI.** Use one background call with a timeout, not one model turn per poll. The CI workflow has 4
jobs (lint and three test jobs, as of 2026-10-08). The Claude workflow adds `claude` check runs, normally
`skipped`, which are not CI:

```bash
sha=$(git rev-parse HEAD)
for i in $(seq 60); do
  out=$(gh api "repos/DD-Ching/jwst-anomaly-research/commits/$sha/check-runs" \
    --jq '[.check_runs[] | select(.name != "claude")] | [length, (map(select(.status != "completed")) | length)] | @tsv')
  [ "$(cut -f1 <<<"$out")" -ge 4 ] && [ "$(cut -f2 <<<"$out")" -eq 0 ] && break
  sleep 30
done
gh api "repos/DD-Ching/jwst-anomaly-research/commits/$sha/check-runs" --jq '.check_runs[] | [.name, .status, .conclusion] | @tsv'
```

Every CI job must conclude `success`; `failure`, `cancelled`, `skipped` or a missing job is not green. If the loop
times out, read `mergeable_state` before waiting again.

**Merging cloud PRs.** CLAUDE.md's merge policy is the gate. From a cloud session the run checks it on GitHub:
- author `DD-Ching` and head repository `DD-Ching/jwst-anomaly-research` (no forks);
- head branch `claude/*`, `batch/*` or `integration/*`;
- no `needs-human` label;
- no guarded file on any page of the PR's files;
- every CI job successful;
- `/code-review` on the final diff, including any conflict resolution;
- `--run-network` tests when I/O code changed;
- no announcement of a result outside the repository.

It then squash-merges with `mcp__github__merge_pull_request` and `expectedHeadSha`. The repository deletes merged
branches itself (`delete_branch_on_merge`). If merging is unavailable or refused, the run labels the PR
`merge-ready` (CI green, review findings fixed, handoff in its CHANGELOG) for you or a local session. Rejected:
a Bash allow rule for the REST merge (#37), because a glob such as `pulls/*/merge` also matches other
`pulls/...` writes (D-028).

**How the PRs reach you.** The run pushes a `claude/<slug>` branch and opens a PR with the `agent` label. Routines
act through your GitHub identity, so these PRs are authored by you. The `agent` label and the session link in the
body tell them apart. Branch protection on `main` requires 0 approvals, so you can merge your own routine PRs. Don't
raise that number, because GitHub doesn't let authors approve their own PRs.

## 4. Desktop scheduled task (local alternative)

This is useful when a run needs local files or the full local disk. In the Desktop app's **Code** tab:
**Routines → New routine → Local**.
- **Instructions:** `/research-cycle`.
- **Folder:** the repository.
- **Worktree toggle on:** each run gets an isolated git worktree.

Click **Run now** once and answer each permission prompt with "always allow", so later runs don't stall. Tasks run
only while the app is open and the computer is awake, and missed runs collapse into a single catch-up run. See
[Desktop scheduled tasks](https://code.claude.com/docs/en/desktop-scheduled-tasks).

## 5. `@claude` in issues and PRs (GitHub Action)

`.github/workflows/claude.yml` runs `anthropics/claude-code-action@v1` when you mention `@claude` in an issue (title
or body), an issue or PR comment, a PR review, or an inline review comment.

Safety:
- Every trigger requires `author_association == 'OWNER'`.
- Only your comments and `claude[bot]`'s are passed to Claude (`include_comments_by_actor`).
- Fork PRs are skipped, so fork code never runs with the token. Someone else's issue or PR body still reaches
  Claude if you mention `@claude` there, so do that only on content you trust.
- Claude pushes a `claude/...` branch, or commits to the PR branch it was invoked on. For issues it posts a link to
  create the PR; it never merges.
- One run per issue or PR at a time. When several `@claude` mentions queue up, for example inline comments in one
  review, GitHub keeps only the newest pending run, and the others show as "cancelled". Put one `@claude` in the
  review summary instead of one per inline comment. The next `/research-cycle` also picks up any owner `@claude`
  comment that `claude[bot]` never answered.
- The `.venv` is built from the default branch. Claude is told to reinstall when the PR changes `pyproject.toml`.

Setup:
```bash
claude.exe setup-token                     # browser login; prints a one-year OAuth token (not saved anywhere)
gh secret set CLAUDE_CODE_OAUTH_TOKEN      # paste the token at the prompt; never commit it
```
Then install the [Claude GitHub App](https://github.com/apps/claude) on the repository. The workflow authenticates
through it with `id-token: write`. Test it with an issue such as `@claude summarize TASKS.md "Now" and propose the
next PR`.

Until the secret exists, the job's first step prints a notice and every later step is skipped, so the workflow is
harmless. The runner installs `.[dev,cloud]` into `.venv` before Claude starts, so Claude can run `python -m pytest`.

The token expires after a year. Set a reminder to rotate it.

Alternative scheduler: an Actions cron (`on: schedule` with `prompt: "/research-cycle"`) also works. It runs only
from the default branch, GitHub disables it after 60 days without repository activity, and it spends Actions
minutes. Prefer the routine.

## 6. Parallel work: when to use `/batch`

The criteria and the ownership rules are in CLAUDE.md "Parallel work": 3 or more independent units, disjoint files,
and a stable interface landed first. When they hold, the cycle fans out instead of suppressing parallelism. How:
- `/batch <instruction>` plans 5-30 units, asks you to approve the plan, then runs one background subagent per unit,
  each in its own git worktree.
- For 3-4 units, or unattended runs where nobody can approve, the cycle uses parallel subagents with
  `isolation: worktree`.
- The batch's PRs share one `batch-<slug>` label and count as one item for the WIP cap. Review them together.

## 7. Your review loop

| Label | Meaning |
|---|---|
| `agent` | Opened by an agent |
| `needs-human` | A decision only you can make: scientific, irreversible, costly or credentials. The PR body says which |
| `candidate` | Candidate report or vetting record |
| `batch-<slug>` | One parallel batch, reviewed together |
| `merge-ready` | A cloud PR that is ready except for the merge itself: merge it |
| `local-wip` | A local session is working on the PR; cloud runs leave it alone |
| `infra`, `science`, `reuse-decision` | Topic |

All these labels already exist in the repository. Agents create new `batch-<slug>` labels with `gh label create`.

1. Review on GitHub with review comments, inline or general.
   - Plain comments are answered by the next `/research-cycle` before any new work (step 2), with a fix commit on the
     PR branch or a reply.
   - For an immediate answer, write `@claude <request>` and the Action handles it. The cycle leaves those comments
     alone unless `claude[bot]` never replied.
2. Agents squash-merge their own PRs once CI is green and `/code-review` is clean, unless the PR carries
   `needs-human` or touches a guarded file (CLAUDE.md "Merge policy and version control"). Those wait for
   you. Other people's PRs always need your approving review. To veto an agent merge after the fact,
   revert the squash commit with a PR.
3. To reject a PR, close it with a one-line reason. The reason is the only record left for the next agent.
4. The queue is TASKS.md. To reprioritize, edit it in a PR, or open an issue yourself. The agents act only on issues
   you authored.

## 8. Safety rails

- **Deny rules** in `.claude/settings.json`, for both Bash and PowerShell:
  - merging around the policy: `gh pr merge --admin` / `--auto` and merge API calls (plain
    `gh pr merge --squash` is allowed for the agents' own PRs);
  - repository settings: branch-protection and ruleset API calls, `gh api` DELETE requests, `gh repo edit`,
    `gh repo delete`;
  - credentials: `gh secret`, `gh auth token`, `gh auth status -t`;
  - force pushes: `--force*`, `-f`, `-uf`, `-qf`, `-vf`, `+refspec`; also `--mirror` and `--all`;
  - pushes that name `main`;
  - remote branch deletion with `--delete` or `-d`;
  - history rewrites: `filter-branch`, `filter-repo`, `reflog expire`;
  - posting local files: `--body-file`, `gh pr/issue -F`, `gh api field=@file` and `--input`, `git diff --no-index`;
  - reading `.env*`, keys and credential files (`~/.ssh`, `~/.aws`, gh's `hosts.yml`, `~/.git-credentials`,
    `~/.claude/.credentials.json`).

  Deny is evaluated before ask and allow, and an allow rule can't override it. These were tested headlessly
  (D-009). The allow list covers the cycle's commands: uv, pytest, ruff, git read, add, commit and switch, pushes to
  `claude/*` and `batch/*`, `gh pr` and `gh issue` reads, comments and creates, and `gh label create`.
- **GitHub MCP tools in cloud sessions** (`mcp__github__*`) are not covered by these Bash rules. Cloud runs use
  them to open and merge their own PRs (D-028). There, the merge policy, branch protection on `main` (required
  checks, enforced for admins) and the repository's merge settings are the boundary.
- **Rules match command text, not intent.** As [the docs warn](https://code.claude.com/docs/en/permissions),
  `git -C . push --force`, `sh -c '...'` and a refspec deletion like `git push origin :branch` slip past them.
  GitHub is the boundary for branches:
  - `main` is protected: PR required, required CI checks, no force pushes, no deletions, enforced for admins.
  - **Recommended:** add a [ruleset](https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/managing-rulesets/about-rulesets)
    on all branches with "Block force pushes" and "Restrict deletions".
  - In cloud runs, the GitHub proxy also rejects branch deletions and tag pushes.
- **Who can merge.** GitHub lets only you and invited collaborators merge. Your agents act as you and merge
  their own PRs under the policy. Required CI checks bind everyone, admins included. Guarded files and
  `needs-human` PRs are left to you by policy, not by GitHub; a separate non-admin bot identity would make
  that a hard boundary (D-009 "Revisit if").
- **Untrusted content:** the skills follow only your instructions (section 1, "Only your words are instructions").
  Text-based deny rules can't prevent every way of leaking a file, so keep secrets off the machines and environments
  where agents run.
- **Action gating:** OWNER-only triggers, filtered comment context, fork PRs skipped, `persist-credentials: false`,
  minimal job permissions. Never enable `show_full_output` or `ACTIONS_STEP_DEBUG`, because the logs of a public
  repository are public.
- **No secrets:** none in the repository, in cloud environment variables or in prompts. A MAST token, if ever
  needed, comes only from env `MAST_API_TOKEN` (CLAUDE.md).
- No hooks yet. They would mainly parse `git push` more robustly, which the recommended ruleset covers server-side
  (D-009).

## 9. Usage and cost

- Routines, `/loop`, Desktop tasks and the Action (through the OAuth token) all draw on your Claude subscription,
  just like interactive use. Check [claude.ai/settings/usage](https://claude.ai/settings/usage).
- Routines also have hourly caps: 100 scheduled runs per hour per account, and 30 "Run now" or API fires per routine
  per hour.
- Cost levers: cadence; the WIP cap (no new PR while 3 await you); `--max-turns` (60 in the Action); the Action's
  60-minute timeout; specific `@claude` requests. Cycles re-read only the state files, not the history.
- GitHub Actions minutes are free for public repositories on standard runners
  ([GitHub billing](https://docs.github.com/en/billing/concepts/product-billing/github-actions)).
- Pause a routine with its on/off switch, rather than deleting it.

## 10. Windows notes

- In Git Bash, call `claude.exe`. The extensionless `claude` shim can be stale. In PowerShell, `claude` works.
- The venv interpreter is `.venv/Scripts/python.exe`; in cloud and Linux it's `.venv/bin/python`.
- The deny rules exist for both the Bash and PowerShell tools, because Claude may use either on Windows.

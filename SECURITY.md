# Security policy

## Supported versions

The project is pre-release (0.0.x). Only the latest commit on `main` receives fixes.

## Reporting a vulnerability

Report privately through GitHub's private vulnerability reporting. In this repository, open the
**Security** tab and choose **Report a vulnerability**, or go directly to
<https://github.com/DD-Ching/jwst-anomaly-research/security/advisories/new>.
Do not open a public issue or pull request for a suspected vulnerability.

Please include the affected file, workflow or commit, the steps to reproduce it, the impact, and a
suggested fix if you have one. The maintainer ([@DD-Ching](https://github.com/DD-Ching)) handles
reports on a best-effort basis. This is a volunteer research project and response times are not
guaranteed. Reporters are credited in the advisory unless they ask not to be.

## Scope

In scope:

- code in `src/`, `scripts/` and `tests/`;
- GitHub Actions workflows and other repository automation, including the AI-agent workflows. For
  example, prompt injection through issues, comments or pull requests that makes an agent leak a
  secret, push code or bypass the owner's review;
- the supply chain: Python dependencies, third-party GitHub Actions and pre-commit hooks;
- secrets or credentials committed by mistake.

Out of scope: vulnerabilities in upstream services or libraries (MAST, astroquery, Astropy, ...),
which should be reported upstream, and scientific disagreements, which belong in a normal issue.

## Data handled

The project processes only public astronomical archive data (JWST/MAST and public catalogs). It does
not collect, store or process personal or user data, runs no hosted service and has no user accounts.

## Secrets policy

- No secrets in the repository, issues, pull requests, logs, notebooks, configs or data manifests.
- A MAST API token, if one is ever needed, is read only from the environment variable
  `MAST_API_TOKEN`. It never goes into code or config files.
- CI and the scheduled network tests run without secrets, using a read-only `GITHUB_TOKEN`. The one
  exception is the failure-report job, which may write issues. Secrets needed by repository automation
  exist only as GitHub Actions secrets, and only the owner manages them.
- If a secret leaks: revoke or rotate it first, then remove it in a new commit and report it privately
  as above. Published git history is not rewritten (project rule), so rotation is the actual fix.

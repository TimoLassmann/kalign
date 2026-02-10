Pre-release check and preparation for version $ARGUMENTS.

Run ALL of the following checks in order. For each check, report PASS or FAIL. Fix issues automatically where possible, and ask before making changes that need judgment.

## 1. Validate version argument

The user must provide a version in `X.Y.Z` format (e.g. `3.4.9`). If `$ARGUMENTS` is empty or not in the right format, stop and ask for the version number.

## 2. Check git state

- Run `git status` — warn about uncommitted changes
- Run `git tag -l "v$ARGUMENTS"` — verify the tag `v$ARGUMENTS` does not already exist
- If the tag exists, warn and ask whether to continue

## 3. Check & fix version numbers

All four files must contain the version `$ARGUMENTS`. Check each one and fix any that don't match:

- **pyproject.toml**: `version = "$ARGUMENTS"`
- **CMakeLists.txt**: `KALIGN_LIBRARY_VERSION_MAJOR`, `KALIGN_LIBRARY_VERSION_MINOR`, `KALIGN_LIBRARY_VERSION_PATCH` must match the three parts of `$ARGUMENTS`
- **build.zig**: `kalignPackageVersion` string and the `KALIGN_PACKAGE_VERSION` in the cflags array must both be `"$ARGUMENTS"`
- **CITATION.cff**: `version:` field must be `$ARGUMENTS`

For any mismatches, fix them directly using the Edit tool and report what was changed.

## 4. Run linting (mirrors CI)

Run these commands in sequence. If `black` or `isort` report formatting issues, auto-fix by running without `--check`:

```
uv run black --check python-kalign/
uv run isort --check-only python-kalign/
uv run flake8 python-kalign --select=E9,F63,F7,F82
```

If black or isort fail the check, run the fix commands:
```
uv run black python-kalign/
uv run isort python-kalign/
```

Then re-run the checks to confirm they pass.

## 5. Run Python tests

```
uv run pytest tests/python/ -v --no-header -q
```

Report the results. Note: there are some known pre-existing test failures in test_parameters.py and test_ecosystem_integration.py — flag any NEW failures beyond those.

## 6. Audit documentation — package name

Search all files for bare `pip install kalign` that should be `pip install kalign-python`. Use grep to find:
- Any `pip install kalign` NOT followed by `-python`
- Specifically check: README files, python-docs/, python-examples/, python-kalign/

If any incorrect references are found, fix them.

## 7. Audit emails

Extract the author email from `pyproject.toml` (the `authors` field). Then search all project files for email addresses and verify they are consistent. Check at minimum:
- `python-kalign/__init__.py` (`__email__`)
- `CITATION.cff`
- `ChangeLog` (most recent entry)

Report any mismatches. Do not hardcode any email — always extract the canonical email from pyproject.toml dynamically.

## 8. Check ChangeLog

Read the `ChangeLog` file and check whether it contains an entry for `version $ARGUMENTS`.

If an entry exists, report PASS.

If NO entry exists:
1. Run `git log $(git describe --tags --abbrev=0)..HEAD --oneline` to get commits since the last tag
2. Read the existing ChangeLog to understand the format (date, author, email, version line, bullet points with changes)
3. Draft a new ChangeLog entry matching the existing format, using today's date and the author name and email from pyproject.toml
4. Present the draft to the user for approval or edits
5. Only after user approval, insert the entry at the top of the ChangeLog file

## 9. Build & verify

Reinstall the package and verify the version:

```
uv pip install -e .
uv run python -c "import kalign; print(kalign.__version__)"
```

Confirm the printed version matches `$ARGUMENTS`. If not, investigate and fix.

## 10. Summary

Print a results table like this:

```
## Release $ARGUMENTS — Pre-release Check Results

| Check                  | Status |
|------------------------|--------|
| Git state              | ...    |
| Version consistency    | ...    |
| Linting (black)        | ...    |
| Linting (isort)        | ...    |
| Linting (flake8)       | ...    |
| Python tests           | ...    |
| Package name audit     | ...    |
| Email consistency      | ...    |
| ChangeLog              | ...    |
| Build & version verify | ...    |
```

Then, if all checks pass (or only known failures remain), ask the user:

> All checks passed. Would you like me to:
> 1. Commit all changes with message "Release vX.Y.Z"
> 2. Create tag `vX.Y.Z`
> 3. Push to remote
>
> Or pick individual steps?

Wait for the user's response before taking any of those actions.

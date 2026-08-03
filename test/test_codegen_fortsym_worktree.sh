#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
storage_root=${FORTFEM_TEST_STORAGE_ROOT:-/mnt/storage}
fixture_root=$(mktemp -d "$storage_root/fortfem-codegen-worktree.XXXXXX")

cleanup() {
    rm -rf "$fixture_root"
}
trap cleanup EXIT

source_repository="$fixture_root/source"
fixture_codegen="$fixture_root/lazy-fortran/fortfem/tools/codegen"
default_worktree="$fixture_root/lazy-fortran/fortsym"
pinned_worktree="$fixture_root/pinned-fortsym"
mismatched_worktree="$fixture_root/mismatched-fortsym"

mkdir -p "$source_repository" "$fixture_codegen"
git -C "$source_repository" init -q
git -C "$source_repository" config user.name "FortFEM test"
git -C "$source_repository" config user.email "fortfem-test@example.invalid"
printf '%s\n' "locked" > "$source_repository/revision.txt"
git -C "$source_repository" add revision.txt
git -C "$source_repository" commit -qm "locked revision"
locked_revision=$(git -C "$source_repository" rev-parse HEAD)
printf '%s\n' "newer" > "$source_repository/revision.txt"
git -C "$source_repository" commit -qam "mismatched revision"
mismatched_revision=$(git -C "$source_repository" rev-parse HEAD)

git -C "$source_repository" worktree add -q --detach \
    "$default_worktree" "$locked_revision"
git -C "$source_repository" worktree add -q --detach \
    "$pinned_worktree" "$locked_revision"
git -C "$source_repository" worktree add -q --detach \
    "$mismatched_worktree" "$mismatched_revision"

cp "$repository_dir/tools/codegen/check_fortsym_revision.sh" "$fixture_codegen/"
cp "$repository_dir/tools/codegen/generate.sh" "$fixture_codegen/"
cp "$repository_dir/tools/codegen/fpm.toml" "$fixture_codegen/"
printf '%s\n' "$locked_revision" > "$fixture_codegen/fortsym.lock"
mkdir -p "$fixture_codegen/app" "$fixture_codegen/src"
checker="$fixture_codegen/check_fortsym_revision.sh"

# The established sibling checkout remains the default.
"$checker"

# An explicit detached worktree at the lock is accepted and remains detached.
FORTSYM_DIR="$pinned_worktree" "$checker"
if git -C "$pinned_worktree" symbolic-ref -q HEAD >/dev/null; then
    echo "pinned FortSym fixture is not detached" >&2
    exit 1
fi

# Generation must build from the override too, rather than merely validating
# it and then resolving the dependency from the default sibling checkout.
fake_bin="$fixture_root/bin"
fake_fo_log="$fixture_root/fake-fo-build-directory"
mkdir -p "$fake_bin"
# The quoted lines are the fake executable's source, not expressions for this
# test process to expand.
# shellcheck disable=SC2016
printf '%s\n' \
    '#!/usr/bin/env bash' \
    'set -euo pipefail' \
    'if [[ "${1:-}" == "build" ]]; then' \
    '    dependency_dir=$(cd ../../../fortsym && pwd -P)' \
    '    [[ "$dependency_dir" == "$FORTSYM_DIR" ]]' \
    '    printf "%s\\n" "$PWD" > "$FORTFEM_FAKE_FO_LOG"' \
    'fi' > "$fake_bin/fo"
chmod +x "$fake_bin/fo"
PATH="$fake_bin:$PATH" \
    FORTSYM_DIR="$pinned_worktree" \
    FORTFEM_FAKE_FO_LOG="$fake_fo_log" \
    FORTFEM_CODEGEN_OUTPUT_DIR="$fixture_root/generated" \
    "$fixture_codegen/generate.sh"
staged_build_directory=$(<"$fake_fo_log")
if [[ "$staged_build_directory" == "$fixture_codegen" ]]; then
    echo "FortSym override did not use an isolated codegen workspace" >&2
    exit 1
fi
if [[ -e "$staged_build_directory" ]]; then
    echo "FortSym override workspace was not cleaned" >&2
    exit 1
fi

run_negative() {
    local label=$1
    local expected=$2
    shift 2
    local output
    if output=$("$@" 2>&1); then
        echo "negative fixture unexpectedly passed: $label" >&2
        return 1
    fi
    grep -Fq "$expected" <<<"$output" || {
        echo "negative fixture $label did not report: $expected" >&2
        printf '%s\n' "$output" >&2
        return 1
    }
}

run_negative mismatched-revision "does not match fortsym.lock" \
    env FORTSYM_DIR="$mismatched_worktree" "$checker"

printf '%s\n' "dirty" >> "$pinned_worktree/revision.txt"
run_negative dirty-worktree "uncommitted changes" \
    env FORTSYM_DIR="$pinned_worktree" "$checker"

echo "FortSym pinned-worktree provenance gate passed"

#!/usr/bin/env just --justfile

# List available recipes
default:
  @just --list

# Get coverage profile
coverage:
  go test -v -coverprofile=provile.cov ./...

# Linting and static checks
lint:
  golangci-lint -c .golangci.yml run

# Run tests
test:
  go test -v ./...
  
branch := `git branch --show-current`
release_changelog := `awk '/## \[Unreleased]/{flag=1; next} /## \[/{flag=0} flag' CHANGELOG.md`

# Bump version in changelog and create commit + tag on current branch
cut-release $NEW_VERSION:
  {{if branch != "release" {error("Releases can only be cut in the releases branch.")} else {""} }}
  just lint
  just test
  sed -i "s/## \[Unreleased]\s*\n*/## [Unreleased]\n\n## [$NEW_VERSION] $(date +'%Y-%m-%d')/g" CHANGELOG.md
  echo "[$NEW_VERSION]: https://github.com/bebop/poly/releases/tag/v$NEW_VERSION" >> CHANGELOG.md
  git add CHANGELOG.md
  git commit -m "Release version v$NEW_VERSION"
  git tag v$NEW_VERSION -m "## Changelog:\n{{release_changelog}}"

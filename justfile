#!/usr/bin/env just --justfile

# List available recipes
default:
  @just --list

# Get coverage profile
coverage:
  go test -v -coverprofile=profile.cov ./...

# Linting and static checks
lint:
  golangci-lint -c .golangci.yml run

# Run tests
test:
  go test -v ./...
  
branch := `git branch --show-current`

# Bump version in changelog and create commit + tag on current branch
cut-release $NEW_VERSION:
  just lint
  just test
  sed -i "s/## \[Unreleased]\s*\n*/## [Unreleased]\n\n## [$NEW_VERSION] $(date +'%Y-%m-%d')/g" CHANGELOG.md
  echo "[$NEW_VERSION]: https://github.com/bebop/poly/releases/tag/v$NEW_VERSION" >> CHANGELOG.md
  git add CHANGELOG.md
  git commit -m "Release version v$NEW_VERSION"
  git tag v$NEW_VERSION

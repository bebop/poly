#!/usr/bin/env just --justfile

default:
  @just --list

coverage:
  go test -v -coverprofile=provile.cov ./...

lint:
  golangci-lint -c .golangci.yml run

test:
  go test -v ./...
  
RELEASE_CHANGELOG := `awk '/## \[Unreleased]/{flag=1; next} /## \[/{flag=0} flag' CHANGELOG.md`

cut-release $NEW_VERSION:
  # Bump version in changelog and create tag for current branch
  sed -i "s/## \[Unreleased]\s*\n*/## [Unreleased]\n\n## [$NEW_VERSION] $(date +'%Y-%m-%d')/g" CHANGELOG.md
  echo "[$NEW_VERSION]: https://github.com/bebop/poly/releases/tag/v$NEW_VERSION" >> CHANGELOG.md
  git add CHANGELOG.md
  git commit -m "Release version v$NEW_VERSION"
  git tag v$NEW_VERSION -m "## Changelog:\n{{RELEASE_CHANGELOG}}"

#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

cd "$(dirname "$0")"/..

go install github.com/golangci/golangci-lint/cmd/golangci-lint@latest
golangci-lint run
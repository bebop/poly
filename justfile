#!/usr/bin/env just --justfile

default:
  @just --list

coverage:
  go test -v -coverprofile=provile.cov ./...

lint:
  golangci-lint -c .golangci.yml run

test:
  just lint
  go test -v ./...
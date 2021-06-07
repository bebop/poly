# See here for image contents: https://github.com/microsoft/vscode-dev-containers/tree/v0.166.1/containers/ubuntu/.devcontainer/base.Dockerfile

# [Choice] Ubuntu version: bionic, focal
ARG VARIANT="focal"
FROM mcr.microsoft.com/vscode/devcontainers/base:0-${VARIANT}


# [Optional] Uncomment this section to install additional OS packages.
RUN export DEBIAN_FRONTEND=noninteractive \
  && sudo apt-get update \
  && sudo apt-get -y install golang

RUN go get -v golang.org/x/tools/gopls
RUN go get -u github.com/lukehoban/go-outline
RUN go get -u github.com/uudashr/gopkgs/cmd/gopkgs
RUN go get github.com/go-delve/delve/cmd/dlv

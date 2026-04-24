#!/usr/bin/env bash
set -euo pipefail

# sracha installer
# Usage:
#   curl -fsSL https://raw.githubusercontent.com/thanhleviet/sracha-rs/main/install.sh | bash
#
# Environment overrides:
#   SRACHA_REPO         GitHub owner/repo to install from (default: thanhleviet/sracha-rs)
#   SRACHA_VERSION      Version tag to install (default: latest release)
#   SRACHA_INSTALL_DIR  Install directory (default: ~/.local/bin)
#   SRACHA_NO_VERIFY    Set to 1 to skip sha256 verification (default: verify when available)

REPO="${SRACHA_REPO:-thanhleviet/sracha-rs}"
BINARY="sracha"
INSTALL_DIR="${SRACHA_INSTALL_DIR:-${HOME}/.local/bin}"

info()  { printf '\033[1;34m==>\033[0m %s\n' "$1"; }
warn()  { printf '\033[1;33mwarning:\033[0m %s\n' "$1" >&2; }
error() { printf '\033[1;31merror:\033[0m %s\n' "$1" >&2; exit 1; }

need() { command -v "$1" >/dev/null 2>&1 || error "required command not found: $1"; }

detect_target() {
    local os arch
    case "$(uname -s)" in
        Linux)  os="linux"  ;;
        Darwin) os="darwin" ;;
        *) error "unsupported OS: $(uname -s). Supported: Linux, macOS." ;;
    esac

    case "$(uname -m)" in
        x86_64|amd64)  arch="x86_64" ;;
        arm64|aarch64) arch="aarch64" ;;
        *) error "unsupported architecture: $(uname -m). Supported: x86_64, aarch64/arm64." ;;
    esac

    case "${os}-${arch}" in
        linux-x86_64)   echo "x86_64-unknown-linux-musl" ;;
        linux-aarch64)  echo "aarch64-unknown-linux-musl" ;;
        darwin-x86_64)  echo "x86_64-apple-darwin" ;;
        darwin-aarch64) echo "aarch64-apple-darwin" ;;
        *) error "no prebuilt binary for ${os}/${arch}" ;;
    esac
}

latest_tag() {
    # Follow the redirect on /releases/latest to avoid the rate-limited API.
    # When a repo has no releases, GitHub redirects to the all-releases page
    # (no /tag/ segment) — return empty in that case so the caller errors cleanly.
    local url
    url=$(curl -fsSLI -o /dev/null -w '%{url_effective}\n' \
        "https://github.com/${REPO}/releases/latest")
    case "$url" in
        *"/releases/tag/"*) printf '%s\n' "${url##*/releases/tag/}" ;;
        *) echo "" ;;
    esac
}

sha256_verify() {
    local archive="$1" sumfile="$2" expected actual
    expected=$(awk '{print $1; exit}' "$sumfile")
    if command -v sha256sum >/dev/null 2>&1; then
        actual=$(sha256sum "$archive" | awk '{print $1}')
    elif command -v shasum >/dev/null 2>&1; then
        actual=$(shasum -a 256 "$archive" | awk '{print $1}')
    else
        warn "no sha256sum/shasum available — skipping checksum verification"
        return 0
    fi
    if [ "$expected" != "$actual" ]; then
        error "checksum mismatch: expected ${expected}, got ${actual}"
    fi
    info "checksum verified (sha256)"
}

update_path_hint() {
    case ":${PATH}:" in
        *":${INSTALL_DIR}:"*) return 0 ;;
    esac

    local shell_rc=""
    case "$(basename "${SHELL:-}")" in
        zsh)  shell_rc="${HOME}/.zshrc"  ;;
        bash) shell_rc="${HOME}/.bashrc" ;;
        *)    [ -f "${HOME}/.profile" ] && shell_rc="${HOME}/.profile" ;;
    esac

    if [ -n "$shell_rc" ] && ! grep -qF "$INSTALL_DIR" "$shell_rc" 2>/dev/null; then
        {
            echo ""
            echo "# Added by sracha installer"
            echo "export PATH=\"${INSTALL_DIR}:\$PATH\""
        } >> "$shell_rc"
        info "added ${INSTALL_DIR} to ${shell_rc}"
        printf '\033[1;33mnotice:\033[0m restart your shell or run:  source %s\n' "$shell_rc"
    else
        warn "${INSTALL_DIR} is not in your PATH"
        echo "  add it manually:"
        echo "    export PATH=\"${INSTALL_DIR}:\$PATH\""
    fi
}

main() {
    need curl
    need tar

    local target tag archive_name url tmp archive sumfile

    target=$(detect_target)
    info "detected target: ${target}"

    if [ -n "${SRACHA_VERSION:-}" ]; then
        tag="${SRACHA_VERSION}"
        info "using specified version: ${tag}"
    else
        info "resolving latest release..."
        tag=$(latest_tag)
        [ -n "$tag" ] || error "could not determine latest version. Set SRACHA_VERSION=vX.Y.Z to install a specific version."
        info "latest version: ${tag}"
    fi

    archive_name="sracha-${target}.tar.gz"
    url="https://github.com/${REPO}/releases/download/${tag}/${archive_name}"

    tmp=$(mktemp -d)
    trap 'rm -rf "${tmp:-}"' EXIT

    archive="${tmp}/${archive_name}"
    sumfile="${archive}.sha256"

    info "downloading ${url}"
    curl -fsSL --retry 3 --retry-delay 2 -o "$archive" "$url" \
        || error "download failed. Verify ${tag} exists at https://github.com/${REPO}/releases"

    if [ "${SRACHA_NO_VERIFY:-0}" != "1" ]; then
        if curl -fsSL --retry 3 --retry-delay 2 -o "$sumfile" "${url}.sha256" 2>/dev/null; then
            ( cd "$tmp" && sha256_verify "$archive_name" "${archive_name}.sha256" )
        else
            warn "no checksum published for ${tag} — skipping verification"
        fi
    fi

    info "extracting..."
    tar -xzf "$archive" -C "$tmp"
    [ -f "${tmp}/${BINARY}" ] || error "binary '${BINARY}' not found in archive"

    mkdir -p "$INSTALL_DIR"
    install -m 0755 "${tmp}/${BINARY}" "${INSTALL_DIR}/${BINARY}"

    info "installed ${BINARY} ${tag} to ${INSTALL_DIR}/${BINARY}"
    echo ""
    "${INSTALL_DIR}/${BINARY}" --version || true
    echo ""

    update_path_hint

    echo ""
    info "try:  ${BINARY} info SRR2584863"
}

main "$@"

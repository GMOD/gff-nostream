# Contributing

## Development

```sh
pnpm install
pnpm test
pnpm build
```

```sh
pnpm version patch  # or minor/major
```

`pnpm version` runs lint, format, types, tests and build, regenerates
CHANGELOG.md with git-cliff, then pushes the tag, which triggers the publish
workflow.

## Publishing

Releases publish automatically via GitHub Actions using npm trusted publishing
(OIDC, no stored token). The workflow requires `--provenance` and
`id-token: write` permissions.

gff-nostream already has trusted publishing set up. To set up a new package:
`npm trust github <pkg> --file publish.yml --repo GMOD/<repo>` (requires
npm >=11.10.0 and 2FA).

Once npm publish succeeds, the `release` job creates the GitHub release for the
tag. `scripts/release-notes.sh` copies the release notes from that version's
CHANGELOG.md section; run it with a version to preview what a release will say.

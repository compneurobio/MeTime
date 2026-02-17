# Release Checklist

> One-line summary: pre-release QA for package and docs website.

## Package

- [ ] `DESCRIPTION` version bumped
- [ ] `NEWS.md` updated
- [ ] `R CMD check` passes

## Documentation

- [ ] README quickstart still valid
- [ ] docs templates updated where needed
- [ ] pkgdown site builds without broken links

## Publish

- [ ] Tag release
- [ ] Publish GitHub release notes
- [ ] Deploy pkgdown site

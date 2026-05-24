# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

**Local development (recommended):**
```bash
docker-compose up
# Access at http://localhost:4000
```

**Local Jekyll directly:**
```bash
bundle install
bundle exec jekyll serve -l -H localhost
# Access at http://localhost:4000
```

**Deploy:** Push to `master` — GitHub Actions auto-builds and publishes to https://DaDa-Chan.github.io.

> Changes to `_config.yml` require restarting the Jekyll server. Other file changes hot-reload automatically.

## Architecture

This is an [Academic Pages](https://academicpages.github.io/) Jekyll site (based on Minimal Mistakes theme) hosted on GitHub Pages.

### Collections (content types)

Content is organized into Jekyll collections. Each collection has a corresponding directory and list page:

| Collection | Directory | List page | Notes |
|---|---|---|---|
| posts (blog) | `_posts/` | `_pages/blog.html` | Filenames must be `YYYY-MM-DD-title.md` |
| plog | `_plog/` | `_pages/plog.html` | Photo journal entries |
| courses | `_courses/` | `_pages/courses.html` | Grouped by `category` front matter field |
| publications | `_publications/` | `_pages/publications.html` | Grouped by `category` (manuscripts/conferences/books) |
| projects | `_projects/` | `_pages/projects.html` | No date prefix required |

### Key configuration files

- **`_config.yml`** — site-wide config: author info, theme, plugins, collection definitions, defaults
- **`_data/navigation.yml`** — top navigation bar menu items
- **`_pages/about.md`** — homepage (permalink: `/`)
- **`_pages/cv.md`** — CV page

### Theming

The active theme is set in `_config.yml` via `site_theme`. Options: `default`, `air`, `sunrise`, `mint`, `dirt`, `contrast`. Each theme has light and dark variants in `_sass/theme/`. CSS variables (e.g. `--global-base-color`, `--global-link-color`) can be overridden in `_includes/head/custom.html` without touching theme files.

Typography variables (`$doc-font-size`, `$type-size-*`, `$sans-serif`) are in `_sass/_themes.scss`.

### MathJax

LaTeX math is supported site-wide via MathJax. Use `$...$` for inline and `$$...$$` for block equations.

### Adding a new collection/section

1. Register in `_config.yml` under `collections:` and add defaults under `defaults:`
2. Create the `_<name>/` directory
3. Create a list page in `_pages/<name>.html` (copy `courses.html` for grouped-by-category layout)
4. Add a nav entry in `_data/navigation.yml`

### Static assets

- Images → `/images/` (referenced as `/images/filename.jpg`)
- Downloadable files (PDFs etc.) → `/files/`
- Profile photo → `/images/profile.png` (controlled by `author.avatar` in `_config.yml`)

"""Check the current public API, compact tutorials, attribution and site links."""
import argparse
import csv
import json
import re
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit


class Page(HTMLParser):
    def __init__(self, path):
        super().__init__()
        self.links, self.words, self.ids = [], [], set()
        self.duplicate_ids, self.images, self.gallery_cards = set(), [], 0
        self.feed(path.read_text())

    def handle_starttag(self, tag, attrs):
        self.links += [v for k, v in attrs if k in ("href", "src") and v]
        for key, value in attrs:
            if key == "id":
                if value in self.ids:
                    self.duplicate_ids.add(value)
                self.ids.add(value)
        if tag == "img":
            self.images.append(dict(attrs))
        if "plot-gallery-card" in dict(attrs).get("class", "").split():
            self.gallery_cards += 1

    def handle_data(self, data):
        self.words.append(data)


def check_site(repo, site):
    repo, site = Path(repo).resolve(), Path(site).resolve()
    errors = []

    def require(condition, message):
        if not condition:
            errors.append(message)

    exports = set(re.findall(r"^export\(([^)]+)\)",
                             (repo / "NAMESPACE").read_text(), re.MULTILINE))
    exports = {name.strip(chr(34) + chr(39)) for name in exports}
    with (repo / "inst/provenance.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    registry = {row["function"]: row for row in rows}
    require(bool(exports), "NAMESPACE has no exports")
    require(exports <= registry.keys(), f"Exports missing provenance: {sorted(exports - registry.keys())}")
    for name in ("inst/NOTICE", "inst/provenance.csv", "CREDITS.md"):
        require((repo / name).is_file() and (repo / name).stat().st_size > 0,
                f"Missing or empty attribution file: {name}")

    aliases = {}
    for path in (repo / "man").glob("*.Rd"):
        for alias in re.findall(r"\\alias\{([^}]+)\}", path.read_text()):
            aliases[alias] = path.stem
    require(exports <= aliases.keys(), f"Exports missing source manuals: {sorted(exports - aliases.keys())}")

    pages = {path.resolve(): Page(path) for path in site.rglob("*.html")}
    for function in sorted(exports & aliases.keys()):
        path = site / "reference" / (aliases[function] + ".html")
        page = pages.get(path)
        require(page is not None, f"Missing reference page: {path.relative_to(site)}")
        if page is not None:
            text = " ".join(" ".join(page.words).split())
            authors = registry.get(function, {}).get("authors", "").split(";")
            require("author" in page.ids and "Author" in text,
                    f"Missing Author section: {path.relative_to(site)}")
            require(any(author.strip() and author.strip() in text for author in authors),
                    f"No credited function author: {path.relative_to(site)}")
            # Full contributor lists may live in the linked CREDITS document.
            require(all(author.strip() in text for author in authors if author.strip()) or
                    any("CREDITS" in urlsplit(link).path for link in page.links),
                    f"Missing full credits or credits link: {path.relative_to(site)}")

    vignette_dir = repo / "vignettes"
    article_sources = [path for path in vignette_dir.rglob("*")
                       if path.suffix.lower() in (".rmd", ".qmd", ".md", ".rnw")]
    tutorials = [path for path in article_sources if path.parent == vignette_dir]
    website_articles = {path.relative_to(vignette_dir).as_posix()
                        for path in article_sources if path.parent != vignette_dir}
    require(len(tutorials) <= 4,
            f"Teaching vignette limit exceeded: {len(tutorials)} > 4")
    approved_articles = {"articles/plot-gallery.Rmd", "articles/ebmf-cnmf-comparison.Rmd",
                         "articles/robustness-checks.Rmd"}
    require(website_articles <= approved_articles,
            f"Unapproved website article: {sorted(website_articles - approved_articles)}")
    articles = {path.stem for path in article_sources}
    require(len(articles) == len(article_sources), "Duplicate article names")
    actual_articles = {path.stem for path in (site / "articles").glob("*.html")} - {"index"}
    require(actual_articles == articles,
            f"Article coverage differs: missing={sorted(articles - actual_articles)}, "
            f"stale={sorted(actual_articles - articles)}")
    require((site / "index.html").is_file(), "Missing site index.html")

    config = (repo / "_pkgdown.yml").read_text()
    match = re.search(r"^url:\s*(\S+)", config, re.MULTILINE)
    canonical = urlsplit(match.group(1).strip(chr(34) + chr(39)) if match else "")
    prefix = canonical.path.rstrip("/") + "/"

    def stale_article(link, source):
        url = urlsplit(link)
        if url.netloc and (url.netloc != canonical.netloc or not url.path.startswith(prefix)):
            return
        match = re.search(r"(?:^|/)articles/([^/?#]+)\.html(?:[?#]|$)", link)
        if match and match.group(1) not in articles | {"index"}:
            errors.append(f"Stale tutorial link: {source}: {link}")

    for path, page in pages.items():
        for link in page.links:
            url = urlsplit(link)
            if url.scheme or url.netloc:
                if not canonical.netloc or url.netloc != canonical.netloc or not url.path.startswith(prefix):
                    continue
            source = path.relative_to(site)
            stale_article(url.path, source)
            target = unquote(url.path)
            if not target:
                dest = path
            elif target.startswith(prefix):
                dest = site / target[len(prefix):]
            elif target.startswith("/"):
                continue
            else:
                dest = path.parent / target
            if dest.is_dir():
                dest = dest / "index.html"
            require(dest.exists(), f"Missing local link: {source}: {link}")
            target_page = pages.get(dest.resolve())
            if url.fragment and target_page is not None:
                require(unquote(url.fragment) in target_page.ids,
                        f"Missing local anchor: {source}: {link}")

    if "plot-gallery" in articles:
        gallery = pages.get(site / "articles/plot-gallery.html")
        if gallery is not None:
            recipes = {name for name in gallery.ids if name.startswith("recipe-")}
            require(gallery.gallery_cards == 17 and len(gallery.images) == 17,
                    "Gallery must have seventeen cards and previews")
            require(len({image.get("src") for image in gallery.images}) == 17,
                    "Gallery must have seventeen distinct preview sources")
            require(len(recipes) == 17, "Gallery must have seventeen recipe targets")
            require(not gallery.duplicate_ids & recipes, "Duplicate gallery recipe ID")
            require(all(image.get("alt", "").strip() for image in gallery.images),
                    "Gallery preview missing alt text")
            linked_recipes = {unquote(urlsplit(link).fragment) for link in gallery.links}
            require(recipes <= linked_recipes, "Unlinked gallery recipe")

    sources = [repo / "README.md", repo / "_pkgdown.yml"]
    sources += list((repo / "R").glob("*.R")) + list((repo / "man").glob("*.Rd"))
    sources += article_sources
    link_pattern = r"(?:https?://[^\s<>]+)?/?articles/[^\s<>]+?\.html"
    for path in sources:
        if path.is_file():
            for link in re.findall(link_pattern, path.read_text()):
                stale_article(link, path.relative_to(repo))
    if errors:
        raise ValueError("\n".join(errors))
    return dict(function_pages=len({aliases[name] for name in exports}),
                exports=len(exports), registry_records=len(rows), articles=len(articles),
                installed_tutorials=len(tutorials), website_articles=len(website_articles),
                author_checks="passed", source_reference_coverage="passed", local_links="passed",
                stale_tutorial_checks="passed", site=str(site))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site", type=Path, default=Path("docs"))
    parser.add_argument("--receipt", type=Path, default=Path("validation/root-site-receipt.json"))
    args = parser.parse_args()
    receipt = check_site(Path(__file__).resolve().parents[2], args.site)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"PASS {receipt['exports']} exports, {receipt['articles']} articles, attribution and local links")


if __name__ == "__main__":
    main()

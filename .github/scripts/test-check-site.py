"""Small independent site fixtures; run with python3 .github/scripts/test-check-site.py."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("check_site", Path(__file__).with_name("check-site.py"))
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


class SiteCheckTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.repo = Path(self.temp.name)
        self.site = self.repo / "site"
        files = {
            "NAMESPACE": "export(foo)\n",
            "inst/provenance.csv": "function,authors\nfoo,Original Author; Other Author\ninternal,Original Author\n",
            "inst/NOTICE": "Retained upstream notice\n",
            "CREDITS.md": "Original Author; Other Author\n",
            "man/foo.Rd": "\\name{foo}\n\\alias{foo}\n",
            "vignettes/beginner.Rmd": "A short tutorial\n",
            "_pkgdown.yml": "url: https://example.org/package/\n",
            "README.md": "[Start](articles/beginner.html)\n",
            "site/index.html": '<a href="articles/beginner.html">Start</a>',
            "site/reference/foo.html": '<h2 id="author">Author</h2><p>Original Author</p><a href="../CREDITS.html">Full credits</a>',
            "site/CREDITS.html": "Original Author; Other Author",
            "site/articles/beginner.html": '<a href="https://example.org/package/reference/foo.html">Function</a>',
        }
        for name, content in files.items():
            path = self.repo / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content)

    def check(self):
        return checker.check_site(self.repo, self.site)

    def test_valid_with_internal_provenance_and_linked_full_credits(self):
        self.assertEqual(self.check()["registry_records"], 2)

    def test_missing_reference_file(self):
        (self.site / "reference/foo.html").unlink()
        with self.assertRaisesRegex(ValueError, "Missing reference page"):
            self.check()

    def test_missing_local_link(self):
        (self.site / "CREDITS.html").unlink()
        with self.assertRaisesRegex(ValueError, "Missing local link"):
            self.check()

    def test_stale_article_even_when_file_still_exists(self):
        (self.site / "articles/obsolete.html").write_text("Old tutorial")
        (self.site / "index.html").write_text('<a href="https://example.org/package/articles/obsolete.html">Old</a>')
        with self.assertRaisesRegex(ValueError, "Stale tutorial link"):
            self.check()

    def test_missing_source_reference(self):
        (self.repo / "man/foo.Rd").unlink()
        with self.assertRaisesRegex(ValueError, "Exports missing source manuals"):
            self.check()

    def test_upstream_articles_are_not_local_tutorials(self):
        for link in ("https://upstream.example/articles/guide.html",
                     "https://example.org/another-package/articles/guide.html"):
            with self.subTest(link=link):
                (self.repo / "vignettes/beginner.Rmd").write_text(f"[Upstream]({link})")
                self.check()

    def test_four_vignette_ceiling(self):
        for i in range(4):
            (self.repo / f"vignettes/extra{i}.Rmd").write_text("Extra")
        with self.assertRaisesRegex(ValueError, "Teaching vignette limit exceeded"):
            self.check()

    def add_gallery(self):
        source = self.repo / "vignettes/articles/plot-gallery.Rmd"
        source.parent.mkdir()
        source.write_text("Website gallery")
        for i in range(12):
            (self.site / f"preview{i}.png").write_bytes(f"image fixture {i}".encode())
        cards = [f'<div class="plot-gallery-card"><a href="#recipe-{i}">'
                 f'<img src="../preview{i}.png" alt="Example {i}"></a></div>'
                 f'<h3 id="recipe-{i}">Recipe {i}</h3>' for i in range(12)]
        self.gallery = self.site / "articles/plot-gallery.html"
        self.gallery.write_text("".join(cards))

    def test_four_tutorials_and_named_website_gallery(self):
        self.add_gallery()
        for i in range(3):
            (self.repo / f"vignettes/tutorial{i}.Rmd").write_text("Tutorial")
            (self.site / f"articles/tutorial{i}.html").write_text("Tutorial")
        receipt = self.check()
        self.assertEqual(receipt["installed_tutorials"], 4)
        self.assertEqual(receipt["website_articles"], 1)
        self.assertEqual(receipt["articles"], 5)

    def test_gallery_does_not_relax_tutorial_ceiling(self):
        self.add_gallery()
        for i in range(4):
            (self.repo / f"vignettes/extra{i}.Rmd").write_text("Extra")
        with self.assertRaisesRegex(ValueError, "Teaching vignette limit exceeded"):
            self.check()

    def test_named_factorization_case_study(self):
        self.add_gallery()
        (self.repo / "vignettes/articles/ebmf-cnmf-comparison.Rmd").write_text("Case study")
        (self.site / "articles/ebmf-cnmf-comparison.html").write_text("Case study")
        receipt = self.check()
        self.assertEqual(receipt["website_articles"], 2)
        self.assertEqual(receipt["articles"], 3)

    def test_unapproved_nested_article(self):
        self.add_gallery()
        (self.repo / "vignettes/articles/extra.Rmd").write_text("Extra")
        with self.assertRaisesRegex(ValueError, "Unapproved website article"):
            self.check()

    def test_missing_gallery_image(self):
        self.add_gallery()
        (self.site / "preview0.png").unlink()
        with self.assertRaisesRegex(ValueError, "Missing local link"):
            self.check()

    def test_missing_same_page_recipe_anchor(self):
        self.add_gallery()
        self.gallery.write_text(self.gallery.read_text().replace(
            'id="recipe-0"', 'id="removed"'))
        with self.assertRaisesRegex(ValueError, "Missing local anchor"):
            self.check()

    def test_duplicate_recipe_anchor(self):
        self.add_gallery()
        with self.gallery.open("a") as handle:
            handle.write('<p id="recipe-0">Duplicate</p>')
        with self.assertRaisesRegex(ValueError, "Duplicate gallery recipe ID"):
            self.check()

    def test_gallery_requires_twelve_previews(self):
        self.add_gallery()
        self.gallery.write_text(self.gallery.read_text().replace(
            '<img src="../preview0.png" alt="Example 0">', ''))
        with self.assertRaisesRegex(ValueError, "twelve cards and previews"):
            self.check()

    def test_gallery_rejects_duplicate_preview_source(self):
        self.add_gallery()
        self.gallery.write_text(self.gallery.read_text().replace('preview1.png', 'preview0.png'))
        with self.assertRaisesRegex(ValueError, "distinct preview sources"):
            self.check()

    def test_gallery_requires_alt_text(self):
        self.add_gallery()
        self.gallery.write_text(self.gallery.read_text().replace('alt="Example 0"', 'alt=""'))
        with self.assertRaisesRegex(ValueError, "missing alt text"):
            self.check()

    def test_relative_and_canonical_anchor_checks(self):
        target = self.site / "articles/beginner.html"
        target.write_text('<h2 id="named-section">Section</h2>')
        index = self.site / "index.html"
        for link in ("articles/beginner.html#named%2Dsection",
                     "https://example.org/package/articles/beginner.html#named-section"):
            with self.subTest(link=link):
                index.write_text(f'<a href="{link}">Section</a>')
                self.check()
                index.write_text(f'<a href="{link}-missing">Section</a>')
                with self.assertRaisesRegex(ValueError, "Missing local anchor"):
                    self.check()


if __name__ == "__main__":
    unittest.main()

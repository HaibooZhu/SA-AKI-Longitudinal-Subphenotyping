#!/usr/bin/env python3
"""Run and archive the PubMed searches used for the JTIM major revision."""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import date
from pathlib import Path


SEARCH_DATE = date(2026, 8, 5)
QUERIES = {
    "sa_aki_trajectory_subphenotype": (
        '("sepsis-associated acute kidney injury"[Title/Abstract] OR '
        '"septic acute kidney injury"[Title/Abstract] OR '
        '(sepsis[Title/Abstract] AND "acute kidney injury"[Title/Abstract])) AND '
        '(trajectory[Title/Abstract] OR trajectories[Title/Abstract] OR '
        'subphenotype*[Title/Abstract] OR phenotype*[Title/Abstract]) AND '
        '("2023/01/01"[Date - Publication] : "3000"[Date - Publication])'
    ),
    "critical_illness_aki_trajectory": (
        '("acute kidney injury"[Title/Abstract]) AND '
        '(trajectory[Title/Abstract] OR trajectories[Title/Abstract] OR '
        'subphenotype*[Title/Abstract]) AND '
        '(critical care[Title/Abstract] OR critically ill[Title/Abstract] OR '
        'intensive care[Title/Abstract]) AND '
        '("2023/01/01"[Date - Publication] : "3000"[Date - Publication])'
    ),
    "furosemide_stress_or_response_aki": (
        '("acute kidney injury"[Title/Abstract]) AND '
        '("furosemide stress test"[Title/Abstract] OR '
        '"furosemide responsiveness"[Title/Abstract] OR '
        '"diuretic responsiveness"[Title/Abstract])'
    ),
}


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo / "results/revision/literature_update",
    )
    return parser.parse_args()


def fetch_json(endpoint: str, params: dict[str, str]) -> dict:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/" + endpoint
    req = urllib.request.Request(
        url + "?" + urllib.parse.urlencode(params),
        headers={"User-Agent": "JTIM-revision-audit/1.0 (literature verification)"},
    )
    with urllib.request.urlopen(req, timeout=60) as response:
        return json.loads(response.read())


def fetch_xml(endpoint: str, params: dict[str, str]) -> ET.Element:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/" + endpoint
    req = urllib.request.Request(
        url + "?" + urllib.parse.urlencode(params),
        headers={"User-Agent": "JTIM-revision-audit/1.0 (literature verification)"},
    )
    with urllib.request.urlopen(req, timeout=60) as response:
        return ET.fromstring(response.read())


def text_or_empty(node: ET.Element | None) -> str:
    return "" if node is None else "".join(node.itertext()).strip()


def parse_article(article: ET.Element) -> dict[str, str]:
    citation = article.find("MedlineCitation")
    pmid = text_or_empty(citation.find("PMID"))
    art = citation.find("Article")
    title = text_or_empty(art.find("ArticleTitle"))
    journal = text_or_empty(art.find("Journal/Title"))
    year = text_or_empty(art.find("Journal/JournalIssue/PubDate/Year"))
    if not year:
        year = text_or_empty(art.find("Journal/JournalIssue/PubDate/MedlineDate"))[:4]
    authors = []
    for author in art.findall("AuthorList/Author"):
        collective = text_or_empty(author.find("CollectiveName"))
        if collective:
            authors.append(collective)
        else:
            surname = text_or_empty(author.find("LastName"))
            initials = text_or_empty(author.find("Initials"))
            if surname:
                authors.append((surname + " " + initials).strip())
    doi = ""
    for item in article.findall("PubmedData/ArticleIdList/ArticleId"):
        if item.attrib.get("IdType") == "doi":
            doi = text_or_empty(item)
            break
    pub_types = "; ".join(text_or_empty(x) for x in art.findall("PublicationTypeList/PublicationType"))
    language = "; ".join(text_or_empty(x) for x in art.findall("Language"))
    return {
        "pmid": pmid,
        "year": year,
        "title": title,
        "authors": "; ".join(authors),
        "journal": journal,
        "doi": doi,
        "publication_types": pub_types,
        "language": language,
        "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
    }


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    archive: dict[str, object] = {
        "search_date": SEARCH_DATE.isoformat(),
        "database": "PubMed via NCBI E-utilities",
        "eligibility": {
            "include": [
                "adult ICU or critical-care AKI/SA-AKI trajectory or subphenotype studies",
                "consensus or methods papers directly relevant to SA-AKI definition or furosemide-response interpretation",
                "English-language peer-reviewed records indexed by the search date",
            ],
            "exclude": [
                "pediatric-only or non-human studies",
                "biomarker-only studies without longitudinal trajectory/subphenotype relevance",
                "conference abstracts, protocols without results, and preprints",
            ],
        },
        "queries": {},
    }
    combined: dict[str, dict[str, str]] = {}
    for name, term in QUERIES.items():
        search = fetch_json(
            "esearch.fcgi",
            {"db": "pubmed", "term": term, "retmode": "json", "retmax": "500", "sort": "pub date"},
        )["esearchresult"]
        ids = search.get("idlist", [])
        records = []
        if ids:
            root = fetch_xml(
                "efetch.fcgi",
                {"db": "pubmed", "id": ",".join(ids), "retmode": "xml"},
            )
            records = [parse_article(x) for x in root.findall("PubmedArticle")]
        for record in records:
            combined[record["pmid"]] = record
        archive["queries"][name] = {
            "term": term,
            "count": int(search["count"]),
            "retrieved": len(records),
            "pmids": ids,
        }
        time.sleep(0.4)

    ordered = sorted(combined.values(), key=lambda x: (x["year"], x["pmid"]), reverse=True)
    archive["unique_records_retrieved"] = len(ordered)
    with (output_dir / "pubmed_search_archive.json").open("w", encoding="utf-8") as handle:
        json.dump(archive, handle, ensure_ascii=False, indent=2)
    with (output_dir / "pubmed_search_records.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(ordered[0]) if ordered else ["pmid"],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(ordered)

    lines = [
        "# Reproducible literature-search record",
        "",
        f"- Search date: {SEARCH_DATE.isoformat()}",
        "- Database: PubMed (NCBI E-utilities)",
        f"- Unique records retrieved: {len(ordered)}",
        "- Screening status: title/abstract relevance screening required before citation insertion",
        "",
        "## Exact queries",
        "",
    ]
    for name, item in archive["queries"].items():
        lines.extend([f"### {name}", "", f"```text\n{item['term']}\n```", "", f"Hits: {item['count']}; retrieved: {item['retrieved']}", ""])
    lines.extend([
        "## Eligibility criteria",
        "",
        "Include adult ICU/critical-care AKI or SA-AKI trajectory/subphenotype studies and directly relevant consensus or furosemide-response methods papers. Exclude pediatric-only, non-human, biomarker-only without trajectory relevance, conference-only, protocol-only, and preprint records.",
        "",
        "The CSV preserves all retrieved records. Citation decisions must be documented separately after relevance screening; retrieval alone does not establish eligibility.",
    ])
    (output_dir / "W7_LITERATURE_SEARCH_PROTOCOL.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()

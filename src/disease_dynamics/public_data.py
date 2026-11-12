"""Cached, polite access to public databases for Track E.

Every response is cached on disk with its URL, parameters, access time and a SHA-256 of
the body, so a rerun is offline and every number can be traced to the query that produced
it. Requests are rate-limited and retried with exponential backoff. Responses are checked
against the fields they are expected to carry, and the module raises instead of guessing.

Cache location: the ``SPATIOTEMPORAL_CACHE`` environment variable, or
``~/.cache/spatiotemporal``. Both are outside OneDrive by default.

Endpoints, as checked on 2026-10-14:

- GTEx Portal API v2 (gtexportal.org/api/v2): gene reference and median expression
- gnomAD GraphQL API (gnomad.broadinstitute.org/api): gene constraint
- GWAS Catalog REST API v2 (www.ebi.ac.uk/gwas/rest/api/v2): associations and metadata
- NCBI E-utilities (eutils.ncbi.nlm.nih.gov): ClinVar record counts
- Human Protein Atlas (www.proteinatlas.org/<ENSG>.json): tissue and cell-type summaries
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import requests

_LAST_REQUEST_TIME = 0.0
# NCBI asks for at most 3 requests per second without an API key; this stays under it
# for every service.
_MIN_INTERVAL_S = 0.4


def cache_directory() -> Path:
    """Return the cache directory, creating it if needed."""
    path = Path(os.environ.get("SPATIOTEMPORAL_CACHE", Path.home() / ".cache" / "spatiotemporal"))
    path.mkdir(parents=True, exist_ok=True)
    return path


@dataclass(frozen=True)
class Provenance:
    """Where a response came from and when."""

    url: str
    accessed_utc: str
    response_sha256: str
    from_cache: bool


def fetch_json(
    url: str,
    params: dict | None = None,
    json_body: dict | None = None,
    timeout_s: float = 30.0,
    retries: int = 4,
) -> tuple[dict, Provenance]:
    """GET (or POST, with ``json_body``) a JSON endpoint through the on-disk cache.

    Parameters
    ----------
    url : str
        Endpoint URL.
    params : dict, optional
        Query parameters.
    json_body : dict, optional
        JSON body; if given, the request is a POST.
    timeout_s : float
        Per-request timeout, seconds.
    retries : int
        Retries after the first attempt, with backoff of 2, 4, 8, ... seconds.

    Returns
    -------
    data : dict
        Parsed JSON.
    provenance : Provenance
        URL, access time and body checksum.
    """
    global _LAST_REQUEST_TIME
    key_text = json.dumps({"url": url, "params": params, "body": json_body}, sort_keys=True)
    cache_file = cache_directory() / f"{hashlib.sha256(key_text.encode()).hexdigest()}.json"
    if cache_file.exists():
        record = json.loads(cache_file.read_text(encoding="utf-8"))
        return record["data"], Provenance(
            record["url"], record["accessed_utc"], record["response_sha256"], True
        )

    last_error: Exception | None = None
    for attempt in range(retries + 1):
        wait = _MIN_INTERVAL_S - (time.monotonic() - _LAST_REQUEST_TIME)
        if wait > 0:
            time.sleep(wait)
        try:
            _LAST_REQUEST_TIME = time.monotonic()
            if json_body is None:
                response = requests.get(url, params=params, timeout=timeout_s)
            else:
                response = requests.post(url, json=json_body, timeout=timeout_s)
            if response.status_code == 429 or response.status_code >= 500:
                raise requests.HTTPError(f"HTTP {response.status_code}", response=response)
            response.raise_for_status()
            body = response.content
            data = response.json()
            break
        except (requests.ConnectionError, requests.Timeout, requests.HTTPError) as error:
            last_error = error
            status = getattr(getattr(error, "response", None), "status_code", None)
            if status is not None and 400 <= status < 500 and status != 429:
                raise  # a client error will not fix itself on retry
            if attempt < retries:
                time.sleep(2.0 ** (attempt + 1))
    else:
        raise ConnectionError(f"{url} failed after {retries + 1} attempts") from last_error

    record = {
        "url": response.url,
        "params": params,
        "body": json_body,
        "accessed_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "response_sha256": hashlib.sha256(body).hexdigest(),
        "data": data,
    }
    # newline="" so Windows does not rewrite line endings in the cached text.
    with open(cache_file, "w", encoding="utf-8", newline="") as handle:
        json.dump(record, handle)
    return data, Provenance(record["url"], record["accessed_utc"], record["response_sha256"], False)


# --- Identifiers and genome builds -------------------------------------------------------

# Every spelling seen in the services used here. An unknown spelling raises rather than
# being guessed.
_BUILD_ALIASES = {
    "GRCh38": "GRCh38",
    "GRCh38/hg38": "GRCh38",
    "GRCh38.p13": "GRCh38",
    "GRCh38.p14": "GRCh38",
    "hg38": "GRCh38",
    "GRCh37": "GRCh37",
    "GRCh37/hg19": "GRCh37",
    "hg19": "GRCh37",
}


def normalise_build(build: str) -> str:
    """Map a genome-build label to GRCh37 or GRCh38, or raise if it is not recognised."""
    if build not in _BUILD_ALIASES:
        raise ValueError(f"unrecognised genome build label {build!r}")
    return _BUILD_ALIASES[build]


def require_same_build(builds: dict[str, str]) -> str:
    """Raise unless every source reports the same genome build; return that build.

    Parameters
    ----------
    builds : dict
        Source name -> build label as the source reports it.
    """
    normalised = {source: normalise_build(label) for source, label in builds.items()}
    if len(set(normalised.values())) != 1:
        raise ValueError(f"genome builds disagree, refusing to join positions: {normalised}")
    return next(iter(normalised.values()))


def stable_ensembl_id(identifier: str) -> str:
    """Strip the version suffix from an Ensembl gene ID (ENSG00000055118.14 -> ENSG00000055118)."""
    if not identifier.startswith("ENSG"):
        raise ValueError(f"not an Ensembl gene ID: {identifier!r}")
    return identifier.split(".")[0]


# --- Service queries ---------------------------------------------------------------------


@dataclass(frozen=True)
class GtexGene:
    gencode_id: str
    symbol: str
    genome_build: str
    gencode_version: str
    provenance: Provenance


def gtex_gene(symbol: str) -> GtexGene:
    """Look up a gene in the GTEx reference; raise unless exactly one match."""
    data, provenance = fetch_json(
        "https://gtexportal.org/api/v2/reference/gene", params={"geneId": symbol}
    )
    matches = [g for g in data.get("data", []) if g.get("geneSymbol") == symbol]
    if len(matches) != 1:
        raise ValueError(f"GTEx returned {len(matches)} genes for {symbol!r}, expected 1")
    gene = matches[0]
    return GtexGene(
        gene["gencodeId"], gene["geneSymbol"], gene["genomeBuild"], gene["gencodeVersion"],
        provenance,
    )


def gtex_median_expression(gencode_id: str, dataset_id: str = "gtex_v8") -> tuple[dict, Provenance]:
    """Median expression per tissue, TPM. Raises if any row reports a different unit."""
    data, provenance = fetch_json(
        "https://gtexportal.org/api/v2/expression/medianGeneExpression",
        params={"gencodeId": gencode_id, "datasetId": dataset_id},
    )
    rows = data.get("data", [])
    if not rows:
        raise ValueError(f"GTEx returned no expression rows for {gencode_id}")
    units = {row["unit"] for row in rows}
    if units != {"TPM"}:
        raise ValueError(f"expected TPM for every tissue, got units {units}")
    if {row["gencodeId"] for row in rows} != {gencode_id}:
        raise ValueError("GTEx rows refer to a different gene ID than requested")
    return {row["tissueSiteDetailId"]: float(row["median"]) for row in rows}, provenance


def gnomad_constraint(symbol: str, reference_genome: str = "GRCh38") -> tuple[dict, Provenance]:
    """Loss-of-function and missense constraint for a gene from gnomAD."""
    query = """
    query Constraint($symbol: String!, $genome: ReferenceGenomeId!) {
      gene(gene_symbol: $symbol, reference_genome: $genome) {
        gene_id
        symbol
        gnomad_constraint { exp_lof obs_lof oe_lof oe_lof_lower oe_lof_upper pli mis_z }
      }
    }"""
    data, provenance = fetch_json(
        "https://gnomad.broadinstitute.org/api",
        json_body={"query": query, "variables": {"symbol": symbol, "genome": reference_genome}},
    )
    if data.get("errors"):
        raise ValueError(f"gnomAD query failed: {data['errors']}")
    gene = data["data"]["gene"]
    if gene is None or gene["gnomad_constraint"] is None:
        raise ValueError(f"gnomAD has no constraint record for {symbol!r}")
    return {"gene_id": gene["gene_id"], "reference_genome": reference_genome,
            **gene["gnomad_constraint"]}, provenance


def clinvar_count(term: str) -> tuple[int, Provenance]:
    """Number of ClinVar records matching an Entrez query. Raises on a mistranslated query."""
    data, provenance = fetch_json(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        params={"db": "clinvar", "term": term, "retmode": "json", "retmax": 0},
    )
    result = data["esearchresult"]
    if "errorlist" in result or "warninglist" in result:
        raise ValueError(f"ClinVar query problem for {term!r}: {result}")
    return int(result["count"]), provenance


def gwas_metadata() -> tuple[dict, Provenance]:
    """GWAS Catalog release metadata, including the genome build of its positions."""
    return fetch_json("https://www.ebi.ac.uk/gwas/rest/api/v2/metadata")


def gwas_associations(mapped_gene: str, page_size: int = 100) -> tuple[list[dict], list[Provenance]]:
    """All GWAS Catalog associations mapped to a gene, following pagination.

    Raises if the number collected differs from the total the API declares.
    """
    associations, provenances, page = [], [], 0
    while True:
        data, provenance = fetch_json(
            "https://www.ebi.ac.uk/gwas/rest/api/v2/associations",
            params={"mapped_gene": mapped_gene, "size": page_size, "page": page},
        )
        provenances.append(provenance)
        for item in data.get("_embedded", {}).get("associations", []):
            associations.append(
                {
                    "rsid": item["snp_allele"][0]["rs_id"] if item.get("snp_allele") else None,
                    "traits": [t["efo_trait"] for t in item.get("efo_traits", [])],
                    "p_value": float(item["p_value"]),
                    "locations": item.get("locations", []),
                }
            )
        page_info = data["page"]
        page += 1
        if page >= page_info["totalPages"]:
            break
    if len(associations) != page_info["totalElements"]:
        raise ValueError(
            f"collected {len(associations)} associations, API declares "
            f"{page_info['totalElements']}"
        )
    return associations, provenances


def hpa_summary(ensembl_id: str) -> tuple[dict, Provenance]:
    """Tissue and single-cell specificity summaries from the Human Protein Atlas."""
    data, provenance = fetch_json(f"https://www.proteinatlas.org/{stable_ensembl_id(ensembl_id)}.json")
    if stable_ensembl_id(data.get("Ensembl", "ENSG")) != stable_ensembl_id(ensembl_id):
        raise ValueError("HPA returned a record for a different gene")
    keys = [
        "Gene",
        "Ensembl",
        "RNA tissue specificity",
        "RNA tissue distribution",
        "RNA tissue specific nTPM",
        "RNA single cell type specificity",
        "RNA single cell type specific nCPM",
        "RNA brain regional specificity",
    ]
    return {key: data.get(key) for key in keys}, provenance

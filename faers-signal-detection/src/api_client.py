"""
api_client.py
-------------
Rate-limited client for the openFDA Drug Adverse Event (FAERS) API.
API docs: https://open.fda.gov/apis/drug/event/
"""

import time
import logging
import requests
from typing import Optional

logger = logging.getLogger(__name__)


class FAERSClient:
    """
    Wrapper around the openFDA /drug/event endpoint.

    Provides methods to query total report counts, drug-specific
    report counts, top adverse events, and background (population)
    adverse event counts — the four values needed to build the
    2x2 contingency table for disproportionality analysis.
    """

    BASE_URL = "https://api.fda.gov/drug/event.json"

    def __init__(self, api_key: Optional[str] = None, delay: float = 0.6):
        """
        Parameters
        ----------
        api_key : str, optional
            openFDA API key (increases rate limit from 40 → 240 req/min).
            Register at https://open.fda.gov/apis/authentication/
        delay : float
            Seconds to wait between requests (default 0.6 s → ~100 req/min).
        """
        self.api_key = api_key
        self.delay = delay
        self.session = requests.Session()

    # ── internal ───────────────────────────────────────────────────────────────

    def _get(self, params: dict) -> Optional[dict]:
        """Make a rate-limited GET to openFDA. Returns None on 404."""
        if self.api_key:
            params["api_key"] = self.api_key
        time.sleep(self.delay)

        try:
            resp = self.session.get(self.BASE_URL, params=params, timeout=30)
            if resp.status_code == 404:
                return None          # no results for this query
            resp.raise_for_status()
            return resp.json()
        except requests.exceptions.RequestException as exc:
            logger.error("openFDA request failed: %s", exc)
            raise

    def _total(self, search: str) -> int:
        """Return total report count for a search string, or 0 if no results."""
        data = self._get({"search": search, "limit": 1})
        if data is None:
            return 0
        return data["meta"]["results"]["total"]

    # ── public ─────────────────────────────────────────────────────────────────

    def get_total_reports(self) -> int:
        """Total adverse event reports in FAERS (background N)."""
        data = self._get({"limit": 1})
        return data["meta"]["results"]["total"]

    def get_drug_total(self, drug_name: str) -> int:
        """Total FAERS reports mentioning *drug_name* (a + b)."""
        return self._total(
            f'patient.drug.openfda.generic_name:"{drug_name}"'
        )

    def get_top_adverse_events(self, drug_name: str, limit: int = 25) -> list:
        """
        Top *limit* adverse events (MedDRA PT) for *drug_name*, each with
        its co-report count with the drug (the 'a' cell of the 2x2 table).
        """
        params = {
            "search": f'patient.drug.openfda.generic_name:"{drug_name}"',
            "count": "patient.reaction.reactionmeddrapt.exact",
            "limit": limit,
        }
        data = self._get(params)
        return data.get("results", []) if data else []

    def get_ae_background_total(self, ae_term: str) -> int:
        """Total FAERS reports mentioning *ae_term* across ALL drugs (a + c)."""
        ae_safe = ae_term.replace('"', '\\"')
        return self._total(
            f'patient.reaction.reactionmeddrapt.exact:"{ae_safe}"'
        )

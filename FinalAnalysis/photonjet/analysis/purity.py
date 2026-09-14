"""Event-leading ABCD region occupancy for PhotonJetTrees_v1."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

import numpy as np
import uproot

from photonjet.analysis.selection import compile_abcd_selection


def _binary_witness(value: Any, *, label: str, path: Path, row: int) -> bool:
    numeric = float(value)
    if not np.isfinite(numeric) or numeric not in {0.0, 1.0}:
        raise ValueError(f"{label} is not binary in {path} row {row}")
    return bool(int(numeric))


def _integer_witness(
    value: Any,
    *,
    label: str,
    path: Path,
    row: int,
    maximum: int,
) -> int:
    if not isinstance(value, (int, np.integer)):
        raise ValueError(f"{label} is not an integer in {path} row {row}")
    numeric = int(value)
    if numeric < 0 or numeric > maximum:
        raise ValueError(f"{label} is outside its range in {path} row {row}")
    return numeric


def purity_counts(
    input_paths: Iterable[Path],
    *,
    non_tight_definition: str,
    isolation_radius: float = 0.4,
) -> dict[str, Any]:
    """Return weighted event-leading A/B/C/D counts.

    Region boundaries are strict: isolated is below the isolated threshold and
    non-isolated is above the non-isolated threshold. Candidates exactly on a
    boundary are in neither region.
    """

    if non_tight_definition not in {"bounded", "complement"}:
        raise ValueError("non_tight_definition must be bounded or complement")
    if isolation_radius not in {0.3, 0.4}:
        raise ValueError("isolation_radius must be 0.3 or 0.4")
    radius = "r03" if isolation_radius == 0.3 else "r04"
    programs = {
        region: compile_abcd_selection(
            region=region,
            non_tight_definition=non_tight_definition,
            isolation_radius=isolation_radius,
        )
        for region in "ABCD"
    }
    branches = [
        "source_file_index",
        "event_id_hi",
        "event_id_lo",
        "photon_encounter_ordinal",
        "photon_et",
        "photon_eta",
        "event_weight",
        "terminal_status",
        "bdt_score",
        "bdt_tight_threshold",
        "bdt_is_tight",
        f"iso_{radius}",
        f"iso_{radius}_threshold",
        f"iso_{radius}_nonisolated_threshold",
        f"iso_{radius}_pass",
    ]
    if non_tight_definition == "bounded":
        branches.extend(
            [
                "bdt_nontight_low_threshold",
                "bdt_nontight_high_threshold",
                "bdt_is_nontight",
            ]
        )
    else:
        branches.append("bdt_is_not_tight")
    counts = np.zeros(4, dtype=float)
    event_counts = np.zeros(4, dtype=int)
    current_key: tuple[int, int, int] | None = None
    current_weight: float | None = None
    current_terminal_status: int | None = None
    current_ordinals: set[int] = set()
    current_leaders: list[tuple[tuple[float, int], float] | None] = [
        None,
        None,
        None,
        None,
    ]
    completed_events: set[tuple[int, int, int]] = set()

    def flush_event() -> None:
        for region_index, leader in enumerate(current_leaders):
            if leader is not None:
                counts[region_index] += leader[1]
                event_counts[region_index] += 1

    for raw_path in input_paths:
        if current_key is not None:
            flush_event()
            completed_events.add(current_key)
            current_key = None
            current_weight = None
            current_terminal_status = None
            current_ordinals = set()
            current_leaders = [None, None, None, None]
        path = Path(raw_path).resolve()
        with uproot.open(path) as root:
            if "photons" not in root:
                raise ValueError(f"{path} does not contain PhotonJetTrees_v1/photons")
            tree = root["photons"]
            missing = sorted(set(branches) - set(tree.keys()))
            if missing:
                raise ValueError(f"{path} photons lacks branches: {', '.join(missing)}")
            arrays = tree.arrays(branches, library="np")
        for index in range(len(arrays["photon_et"])):
            terminal_status = arrays["terminal_status"][index]
            event_accepted = programs["A"].accepts_event(terminal_status)
            photon_et = float(arrays["photon_et"][index])
            photon_eta = float(arrays["photon_eta"][index])
            encounter_ordinal = _integer_witness(
                arrays["photon_encounter_ordinal"][index],
                label="photon encounter ordinal",
                path=path,
                row=index,
                maximum=2**31 - 1,
            )
            event_weight = float(arrays["event_weight"][index])
            if (
                not np.all(np.isfinite([photon_et, photon_eta, event_weight]))
            ):
                raise ValueError(
                    f"non-finite or invalid photon identity in {path} row {index}"
                )
            source_file_index = _integer_witness(
                arrays["source_file_index"][index],
                label="source file index",
                path=path,
                row=index,
                maximum=2**31 - 1,
            )
            event_id_hi = _integer_witness(
                arrays["event_id_hi"][index],
                label="event identity high word",
                path=path,
                row=index,
                maximum=2**64 - 1,
            )
            event_id_lo = _integer_witness(
                arrays["event_id_lo"][index],
                label="event identity low word",
                path=path,
                row=index,
                maximum=2**64 - 1,
            )
            if event_id_hi == 0 and event_id_lo == 0:
                raise ValueError(f"event identity is null in {path} row {index}")
            key = (
                source_file_index,
                event_id_hi,
                event_id_lo,
            )
            if key != current_key:
                if current_key is not None:
                    flush_event()
                    completed_events.add(current_key)
                if key in completed_events:
                    raise ValueError("photon rows for one event are not contiguous")
                current_key = key
                current_weight = event_weight
                current_terminal_status = int(terminal_status)
                current_ordinals = set()
                current_leaders = [None, None, None, None]
            if current_weight != event_weight:
                raise ValueError("photon rows in one event carry different event weights")
            if current_terminal_status != terminal_status:
                raise ValueError("photon rows in one event carry different terminal status")
            if encounter_ordinal in current_ordinals:
                raise ValueError("photon encounter ordinal is duplicated within an event")
            current_ordinals.add(encounter_ordinal)
            if not event_accepted:
                continue
            candidate = {
                "photon_et": photon_et,
                "photon_eta": photon_eta,
                "photon_encounter_ordinal": encounter_ordinal,
            }
            if not programs["A"].accepts("photon", candidate):
                continue

            score = float(arrays["bdt_score"][index])
            tight_threshold = float(arrays["bdt_tight_threshold"][index])
            isolation = float(arrays[f"iso_{radius}"][index])
            isolation_threshold = float(arrays[f"iso_{radius}_threshold"][index])
            nonisolated_threshold = float(
                arrays[f"iso_{radius}_nonisolated_threshold"][index]
            )
            numeric = [
                score,
                tight_threshold,
                isolation,
                isolation_threshold,
                nonisolated_threshold,
            ]
            non_tight_low: float | None = None
            non_tight_high: float | None = None
            if non_tight_definition == "bounded":
                non_tight_low = float(arrays["bdt_nontight_low_threshold"][index])
                non_tight_high = float(arrays["bdt_nontight_high_threshold"][index])
                numeric.extend([non_tight_low, non_tight_high])
            if (
                np.any(~np.isfinite(np.asarray(numeric, dtype=float)))
                or isolation_threshold > nonisolated_threshold
                or (
                    non_tight_definition == "bounded"
                    and not float(non_tight_low)
                    < float(non_tight_high)
                    <= tight_threshold
                )
            ):
                raise ValueError(
                    f"non-finite or invalid photon witnesses in {path} row {index}"
                )
            tight = score > tight_threshold
            isolated = isolation < isolation_threshold
            if _binary_witness(
                arrays["bdt_is_tight"][index],
                label="tight BDT witness",
                path=path,
                row=index,
            ) != tight:
                raise ValueError(
                    f"tight BDT witness differs from its score in {path} row {index}"
                )
            if non_tight_definition == "bounded":
                bounded = float(non_tight_low) < score < float(non_tight_high)
                if _binary_witness(
                    arrays["bdt_is_nontight"][index],
                    label="bounded BDT witness",
                    path=path,
                    row=index,
                ) != bounded:
                    raise ValueError(
                        f"bounded BDT witness differs from its score in {path} row {index}"
                    )
            elif _binary_witness(
                arrays["bdt_is_not_tight"][index],
                label="complement BDT witness",
                path=path,
                row=index,
            ) != (not tight):
                raise ValueError(
                    f"complement BDT witness differs from its score in {path} row {index}"
                )
            if _binary_witness(
                arrays[f"iso_{radius}_pass"][index],
                label=f"isolation {radius} witness",
                path=path,
                row=index,
            ) != isolated:
                raise ValueError(
                    f"isolation witness differs from its value in {path} row {index}"
                )
            candidate.update(
                {
                    "photon_bdt_score": score,
                    "photon_bdt_tight_threshold": tight_threshold,
                    f"photon_iso_{radius}": isolation,
                    f"photon_iso_{radius}_threshold": isolation_threshold,
                    f"photon_iso_{radius}_nonisolated_threshold": (
                        nonisolated_threshold
                    ),
                }
            )
            if non_tight_definition == "bounded":
                candidate["photon_bdt_nontight_low_threshold"] = float(non_tight_low)
                candidate["photon_bdt_nontight_high_threshold"] = float(non_tight_high)
            memberships = [
                region_index
                for region_index, region in enumerate("ABCD")
                if programs[region].accepts("photon", candidate)
                and programs[region].accepts("photon_class", candidate)
            ]
            if len(memberships) > 1:
                raise ValueError("one photon satisfies multiple exclusive ABCD regions")
            if memberships:
                rank = programs["A"].leader_rank(candidate)
                region_index = memberships[0]
                previous = current_leaders[region_index]
                if previous is None or rank < previous[0]:
                    current_leaders[region_index] = (rank, event_weight)
    if current_key is not None:
        flush_event()
    return {
        "schema": "PhotonJetABCDCountsV1",
        "non_tight_definition": non_tight_definition,
        "isolation_radius": isolation_radius,
        "selection_program_sha256": {
            region: program.sha256 for region, program in programs.items()
        },
        "weighted_counts": dict(zip("ABCD", counts.tolist())),
        "event_counts": dict(zip("ABCD", event_counts.tolist())),
    }


def raw_abcd_purity(counts: dict[str, float]) -> float:
    """Compute the uncorrected ABCD estimate for diagnostics.

    Physics results must apply the declared leakage and closure treatment
    before using this quantity.
    """

    a, b, c, d = (float(counts[key]) for key in "ABCD")
    if a <= 0 or d <= 0:
        raise ValueError("raw ABCD purity requires positive A and D")
    return 1.0 - (b * c) / (a * d)

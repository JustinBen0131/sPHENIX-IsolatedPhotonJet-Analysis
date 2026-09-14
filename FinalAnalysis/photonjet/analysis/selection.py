"""Executable selection contracts shared by reducers and plot annotations.

The predicate objects in this module are deliberately small. Each one owns
the field, transformation, comparison, and value used by the event loop. The
same immutable objects are serialized into the histogram receipt and later
formatted by the plotting layer. A cut therefore cannot change without also
changing the contract from which its visible label is compiled.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math
from numbers import Integral
from typing import Any, Mapping, Sequence

from photonjet.provenance import canonical_json


CANONICAL_LEADER_PHOTON_ET_MIN = 15.0
CANONICAL_LEADER_PHOTON_ET_MAX = 35.0
CANONICAL_LEADER_PHOTON_ABS_ETA_MAX = 0.7


@dataclass(frozen=True)
class Predicate:
    """One executable numeric requirement in a named analysis stage."""

    key: str
    stage: str
    field: str
    operator: str
    value: float | None = None
    lower: float | None = None
    upper: float | None = None
    tolerance: float | None = None
    reference: str | None = None
    lower_reference: str | None = None
    upper_reference: str | None = None
    transform: str = "identity"
    unit: str | None = None

    def __post_init__(self) -> None:
        if self.transform not in {"identity", "abs"}:
            raise ValueError(f"unsupported transform: {self.transform}")
        parameter_fields = {
            "value", "lower", "upper", "tolerance", "reference",
            "lower_reference", "upper_reference",
        }
        consumed_parameters = {
            "eq": {"value"},
            "gt": {"value"},
            "lt": {"value"},
            "closed_open": {"lower", "upper"},
            "near": {"value", "tolerance"},
            "gt_field": {"reference"},
            "lt_field": {"reference"},
            "open_between_fields": {"lower_reference", "upper_reference"},
            "not_gt_field": {"reference"},
        }
        if self.operator not in consumed_parameters:
            raise ValueError(f"unsupported operator: {self.operator}")
        supplied = {
            name for name in parameter_fields if getattr(self, name) is not None
        }
        unexpected = sorted(supplied - consumed_parameters[self.operator])
        if unexpected:
            raise ValueError(
                f"{self.operator} has non-executable parameters: {unexpected}"
            )
        if self.operator in {"eq", "gt", "lt", "near"} and self.value is None:
            raise ValueError(f"{self.operator} requires value")
        if self.operator == "closed_open":
            if self.lower is None or self.upper is None or not self.lower < self.upper:
                raise ValueError("closed_open requires an increasing interval")
        if self.operator == "near" and (self.tolerance is None or self.tolerance <= 0):
            raise ValueError("near requires a positive tolerance")
        if self.operator in {"gt_field", "lt_field", "not_gt_field"} and not self.reference:
            raise ValueError(f"{self.operator} requires a reference field")
        if self.operator == "open_between_fields" and not (
            self.lower_reference and self.upper_reference
        ):
            raise ValueError("open_between_fields requires lower and upper reference fields")
        for number in (self.value, self.lower, self.upper, self.tolerance):
            if number is not None and not math.isfinite(number):
                raise ValueError(f"non-finite predicate value for {self.key}")

    def evaluate(self, record: Mapping[str, float]) -> bool:
        if self.field not in record:
            raise KeyError(f"selection field is absent: {self.field}")
        observed = float(record[self.field])
        if not math.isfinite(observed):
            return False
        if self.transform == "abs":
            observed = abs(observed)
        if self.operator == "eq":
            return observed == float(self.value)
        if self.operator == "gt":
            return observed > float(self.value)
        if self.operator == "lt":
            return observed < float(self.value)
        if self.operator == "closed_open":
            return float(self.lower) <= observed < float(self.upper)
        if self.operator == "near":
            return abs(observed - float(self.value)) < float(self.tolerance)
        if self.operator in {"gt_field", "lt_field", "not_gt_field"}:
            if self.reference not in record:
                raise KeyError(f"selection reference field is absent: {self.reference}")
            reference = float(record[self.reference])
            if not math.isfinite(reference):
                return False
            if self.operator == "gt_field":
                return observed > reference
            if self.operator == "lt_field":
                return observed < reference
            return not observed > reference
        if self.operator == "open_between_fields":
            if self.lower_reference not in record or self.upper_reference not in record:
                raise KeyError("selection interval reference field is absent")
            lower = float(record[self.lower_reference])
            upper = float(record[self.upper_reference])
            if not math.isfinite(lower) or not math.isfinite(upper) or not lower < upper:
                return False
            return lower < observed < upper
        raise AssertionError(self.operator)

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {
            "key": self.key,
            "stage": self.stage,
            "field": self.field,
            "operator": self.operator,
            "transform": self.transform,
        }
        for name in (
            "value", "lower", "upper", "tolerance", "reference",
            "lower_reference", "upper_reference", "unit",
        ):
            value = getattr(self, name)
            if value is not None:
                result[name] = value
        return result

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "Predicate":
        allowed = {
            "key", "stage", "field", "operator", "value", "lower", "upper",
            "tolerance", "reference", "lower_reference", "upper_reference",
            "transform", "unit",
        }
        extra = sorted(set(value) - allowed)
        if extra:
            raise ValueError(f"unknown predicate fields: {extra}")
        return cls(**dict(value))


@dataclass(frozen=True)
class SelectionFact:
    """A categorical choice that affects which rows enter the payload."""

    key: str
    value: str | float

    def to_dict(self) -> dict[str, Any]:
        return {"key": self.key, "value": self.value}

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "SelectionFact":
        if set(value) != {"key", "value"}:
            raise ValueError("selection fact must contain exactly key and value")
        return cls(key=str(value["key"]), value=value["value"])


def _canonical_event_predicates() -> tuple[Predicate, ...]:
    """Retained diagnostic rows are not accepted physics events.

    The producer owns the rejection reason. Offline consumers use only its
    recorded success state, without reimplementing detector-specific gates.
    """

    return (
        Predicate(
            key="event.accepted",
            stage="event",
            field="terminal_status",
            operator="eq",
            value=0.0,
        ),
    )


def _canonical_leader_photon_predicates() -> tuple[Predicate, ...]:
    """Eligibility encoded by the recorded V1 event-leader branches."""

    return (
        Predicate(
            key="photon.et",
            stage="photon",
            field="photon_et",
            operator="closed_open",
            lower=CANONICAL_LEADER_PHOTON_ET_MIN,
            upper=CANONICAL_LEADER_PHOTON_ET_MAX,
            unit="GeV",
        ),
        Predicate(
            key="photon.abs_eta",
            stage="photon",
            field="photon_eta",
            operator="lt",
            transform="abs",
            value=CANONICAL_LEADER_PHOTON_ABS_ETA_MAX,
        ),
    )


@dataclass(frozen=True)
class SelectionProgram:
    """The complete, executable lineage for one recoil payload."""

    predicates: tuple[Predicate, ...]
    facts: tuple[SelectionFact, ...]
    leader_branch: str | None

    def __post_init__(self) -> None:
        predicate_keys = [predicate.key for predicate in self.predicates]
        if len(predicate_keys) != len(set(predicate_keys)):
            raise ValueError("duplicate selection predicate key")
        event_predicates = tuple(
            predicate for predicate in self.predicates if predicate.stage == "event"
        )
        if event_predicates != _canonical_event_predicates():
            raise ValueError("selection requires the canonical accepted-event predicate")
        fact_keys = [fact.key for fact in self.facts]
        if len(fact_keys) != len(set(fact_keys)):
            raise ValueError("duplicate selection fact key")
        scopes = [
            fact.value
            for fact in self.facts
            if fact.key == "photon.candidate_scope"
        ]
        if len(scopes) != 1 or scopes[0] not in {"inclusive_pairs", "event_leading"}:
            raise ValueError("selection program requires one supported candidate scope")
        if scopes[0] == "inclusive_pairs":
            if self.facts != (
                SelectionFact("photon.candidate_scope", "inclusive_pairs"),
            ):
                raise ValueError("inclusive selection has non-canonical facts")
            if self.leader_branch is not None:
                raise ValueError("inclusive selection cannot claim a recorded leader branch")
            if any(
                predicate.stage == "photon_class"
                for predicate in self.predicates
            ):
                raise ValueError("inclusive selection cannot carry ABCD predicates")
            return

        fact_values = {fact.key: fact.value for fact in self.facts}
        try:
            region = str(fact_values["photon.abcd_region"])
            non_tight_definition = str(
                fact_values.get("photon.non_tight_definition", "bounded")
            )
            isolation_radius = float(fact_values["photon.isolation_radius"])
            expected_predicates, expected_facts, expected_leader = (
                _canonical_abcd_components(
                    region=region,
                    non_tight_definition=non_tight_definition,
                    isolation_radius=isolation_radius,
                )
            )
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("event-leading selection facts are not canonical") from error
        photon_class_predicates = tuple(
            predicate
            for predicate in self.predicates
            if predicate.stage == "photon_class"
        )
        photon_predicates = tuple(
            predicate
            for predicate in self.predicates
            if predicate.stage == "photon"
        )
        program_leader = (
            expected_leader
            if photon_predicates == _canonical_leader_photon_predicates()
            else None
        )
        if (
            photon_class_predicates != expected_predicates
            or self.facts != expected_facts
            or self.leader_branch != program_leader
        ):
            raise ValueError(
                "event-leading facts, predicates, and leader branch disagree"
            )

    def accepts(self, stage: str, record: Mapping[str, float]) -> bool:
        return all(
            predicate.evaluate(record)
            for predicate in self.predicates
            if predicate.stage == stage
        )

    def accepts_event(self, terminal_status: int) -> bool:
        """Execute event acceptance without coercing corrupt status values.

        In particular, converting a fractional value to int could turn a
        malformed rejection state into zero and admit it as a physics event.
        """

        if (
            isinstance(terminal_status, bool)
            or not isinstance(terminal_status, Integral)
            or not 0 <= terminal_status <= 2**31 - 1
        ):
            raise ValueError("terminal_status must be a non-negative int32")
        return self.accepts("event", {"terminal_status": terminal_status})

    def leader_rank(self, record: Mapping[str, float]) -> tuple[float, int]:
        """Return the one canonical event-leader ordering key."""

        scope = next(
            fact.value
            for fact in self.facts
            if fact.key == "photon.candidate_scope"
        )
        if scope != "event_leading":
            raise ValueError("inclusive selection has no event-leading candidate")
        if "photon_encounter_ordinal" not in record:
            raise KeyError("selection field is absent: photon_encounter_ordinal")
        raw_ordinal = float(record["photon_encounter_ordinal"])
        if (
            not math.isfinite(raw_ordinal)
            or not raw_ordinal.is_integer()
            or raw_ordinal < 0
        ):
            raise ValueError("photon encounter ordinal is not a non-negative integer")
        if "photon_et" not in record:
            raise KeyError("selection field is absent: photon_et")
        photon_et = float(record["photon_et"])
        if not math.isfinite(photon_et):
            raise ValueError("eligible photon ET is not finite")
        return (-photon_et, int(raw_ordinal))

    def choose_leader(self, records: Sequence[Mapping[str, float]]) -> int:
        """Execute the serialized ABCD predicates and deterministic leader rule."""

        scope = next(
            fact.value
            for fact in self.facts
            if fact.key == "photon.candidate_scope"
        )
        if scope != "event_leading":
            raise ValueError("inclusive selection has no event-leading candidate")
        ranks = [self.leader_rank(record) for record in records]
        ordinals = [rank[1] for rank in ranks]
        if len(ordinals) != len(set(ordinals)):
            raise ValueError("photon encounter ordinal is duplicated within an event")
        eligibility_stages = (
            ("photon", "photon_class")
            if any(predicate.stage == "photon" for predicate in self.predicates)
            else ("photon_class",)
        )
        eligible = [
            index
            for index, record in enumerate(records)
            if all(self.accepts(stage, record) for stage in eligibility_stages)
        ]
        if not eligible:
            return -1
        return min(
            eligible,
            key=lambda index: (*ranks[index], index),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": "PhotonJetSelectionProgramV1",
            "predicates": [predicate.to_dict() for predicate in self.predicates],
            "facts": [fact.to_dict() for fact in self.facts],
            "leader_branch": self.leader_branch,
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(canonical_json(self.to_dict()).encode("utf-8")).hexdigest()

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "SelectionProgram":
        if set(value) != {"schema", "predicates", "facts", "leader_branch"}:
            raise ValueError("selection program fields differ from PhotonJetSelectionProgramV1")
        if value["schema"] != "PhotonJetSelectionProgramV1":
            raise ValueError("unsupported selection program schema")
        predicates = value["predicates"]
        facts = value["facts"]
        if not isinstance(predicates, list) or not isinstance(facts, list):
            raise ValueError("selection predicates and facts must be lists")
        leader = value["leader_branch"]
        if leader is not None and not isinstance(leader, str):
            raise ValueError("leader_branch must be a string or null")
        return cls(
            predicates=tuple(Predicate.from_dict(item) for item in predicates),
            facts=tuple(SelectionFact.from_dict(item) for item in facts),
            leader_branch=leader,
        )


def _canonical_abcd_components(
    *,
    region: str,
    non_tight_definition: str,
    isolation_radius: float,
) -> tuple[tuple[Predicate, ...], tuple[SelectionFact, ...], str | None]:
    """Build the one canonical ABCD semantic signature without recursion."""

    if region not in {"A", "B", "C", "D"}:
        raise ValueError(f"unsupported ABCD region: {region}")
    if non_tight_definition not in {"bounded", "complement"}:
        raise ValueError(
            f"unsupported non-tight definition: {non_tight_definition}"
        )
    if isolation_radius not in {0.3, 0.4}:
        raise ValueError("isolation radius must be 0.3 or 0.4")

    identity = "tight" if region in {"A", "B"} else "non_tight"
    isolation = "isolated" if region in {"A", "C"} else "nonisolated"
    radius_key = "r03" if isolation_radius == 0.3 else "r04"
    facts: list[SelectionFact] = [
        SelectionFact("photon.candidate_scope", "event_leading"),
        SelectionFact("photon.abcd_region", region),
        SelectionFact("photon.id_state", identity),
        SelectionFact("photon.isolation_state", isolation),
        SelectionFact("photon.isolation_radius", isolation_radius),
        SelectionFact("photon.score_source", "recorded_bdt_score"),
        SelectionFact(
            "photon.bdt_threshold_source",
            "per_candidate_calibrated_threshold",
        ),
        SelectionFact(
            "photon.isolation_threshold_source",
            "per_candidate_isolation_witness",
        ),
        SelectionFact(
            "photon.leader_rule",
            "highest_et_then_encounter_ordinal",
        ),
    ]
    if identity == "non_tight":
        facts.append(
            SelectionFact("photon.non_tight_definition", non_tight_definition)
        )

    if identity == "tight":
        classification = Predicate(
            key="photon.bdt_class",
            stage="photon_class",
            field="photon_bdt_score",
            operator="gt_field",
            reference="photon_bdt_tight_threshold",
        )
    elif non_tight_definition == "bounded":
        classification = Predicate(
            key="photon.bdt_class",
            stage="photon_class",
            field="photon_bdt_score",
            operator="open_between_fields",
            lower_reference="photon_bdt_nontight_low_threshold",
            upper_reference="photon_bdt_nontight_high_threshold",
        )
    else:
        classification = Predicate(
            key="photon.bdt_class",
            stage="photon_class",
            field="photon_bdt_score",
            operator="not_gt_field",
            reference="photon_bdt_tight_threshold",
        )
    isolation_predicate = Predicate(
        key="photon.isolation_class",
        stage="photon_class",
        field=f"photon_iso_{radius_key}",
        operator="lt_field" if isolation == "isolated" else "gt_field",
        reference=(
            f"photon_iso_{radius_key}_threshold"
            if isolation == "isolated"
            else f"photon_iso_{radius_key}_nonisolated_threshold"
        ),
        unit="GeV",
    )
    suffix = (
        "_complement"
        if region in {"C", "D"} and non_tight_definition == "complement"
        else ""
    )
    leader_branch = (
        f"leader_{region}_r04{suffix}_index"
        if isolation_radius == 0.4
        else None
    )
    return (
        (classification, isolation_predicate),
        tuple(facts),
        leader_branch,
    )


def compile_abcd_selection(
    *,
    region: str,
    non_tight_definition: str,
    isolation_radius: float,
) -> SelectionProgram:
    """Compile the sole executable ABCD classification and leader program."""

    predicates, facts, leader_branch = _canonical_abcd_components(
        region=region,
        non_tight_definition=non_tight_definition,
        isolation_radius=isolation_radius,
    )
    return SelectionProgram(
        predicates=(
            *_canonical_event_predicates(),
            *_canonical_leader_photon_predicates(),
            *predicates,
        ),
        facts=facts,
        leader_branch=leader_branch,
    )


def compile_recoil_selection(selection: Any) -> SelectionProgram:
    """Compile a RecoilSelection-like object into the sole cut program."""

    region = str(selection.region)
    if region not in {"inclusive", "A", "B", "C", "D"}:
        raise ValueError(f"unsupported photon region: {region}")
    non_tight = str(selection.non_tight_definition)
    if non_tight not in {"bounded", "complement"}:
        raise ValueError(f"unsupported non-tight definition: {non_tight}")

    predicates = (
        *_canonical_event_predicates(),
        Predicate(
            key="photon.et",
            stage="photon",
            field="photon_et",
            operator="closed_open",
            lower=float(selection.photon_et_min),
            upper=float(selection.photon_et_max),
            unit="GeV",
        ),
        Predicate(
            key="photon.abs_eta",
            stage="photon",
            field="photon_eta",
            operator="lt",
            transform="abs",
            value=float(selection.photon_abs_eta_max),
        ),
        Predicate(
            key="jet.pt",
            stage="recoil",
            field="jet_pt",
            operator="gt",
            value=float(selection.jet_pt_min),
            unit="GeV",
        ),
        Predicate(
            key="jet.abs_eta",
            stage="recoil",
            field="jet_eta",
            operator="lt",
            transform="abs",
            value=float(selection.jet_abs_eta_max),
        ),
        Predicate(
            key="jet.radius",
            stage="recoil",
            field="jet_radius",
            operator="near",
            value=float(selection.jet_radius),
            tolerance=1.0e-9,
        ),
        Predicate(
            key="recoil.delta_phi",
            stage="recoil",
            field="delta_phi",
            operator="gt",
            value=float(selection.delta_phi_min),
            unit="rad",
        ),
    )

    if region == "inclusive":
        return SelectionProgram(
            predicates=predicates,
            facts=(SelectionFact("photon.candidate_scope", "inclusive_pairs"),),
            leader_branch=None,
        )

    classification_predicates, facts, canonical_leader_branch = (
        _canonical_abcd_components(
            region=region,
            non_tight_definition=non_tight,
            isolation_radius=float(getattr(selection, "isolation_radius", 0.4)),
        )
    )
    photon_predicates = tuple(
        predicate for predicate in predicates if predicate.stage == "photon"
    )
    leader_branch = (
        canonical_leader_branch
        if photon_predicates == _canonical_leader_photon_predicates()
        else None
    )
    return SelectionProgram(
        predicates=(*predicates, *classification_predicates),
        facts=facts,
        leader_branch=leader_branch,
    )

"""
Core sequence prediction from precursor peptides.
License: MIT
Author: Magnus Ohle

This module provides rule-based prediction of lasso peptide core sequences
from full precursor sequences, using configurable rules and patterns
learned from training data.
"""

from typing import Dict, List, Optional, Tuple
from pathlib import Path

from lasso_workbench.schemas.prediction import PredictionResult, PrecursorPrediction, Orientation
from lasso_workbench.schemas.config import RuleSet, RuleParameter
from lasso_workbench.config.rules_io import (
    load_active_ruleset as _load_active_ruleset,
    DEFAULT_RULES_PATH,
)
from Bio.SeqUtils import ProtParamData
from Bio.SeqUtils.ProtParam import ProteinAnalysis

_KYTE_DOOLITTLE = ProtParamData.kd


def load_active_ruleset() -> RuleSet:
    """Load user rules if present, otherwise defaults."""
    return _load_active_ruleset()


def load_default_ruleset() -> RuleSet:
    """Load the default ruleset only."""
    return _load_active_ruleset(Path("__missing__"), DEFAULT_RULES_PATH)


class RuleEngine:
    """
    Pure logic engine for core prediction.
    
    Initialized with a RuleSet schema, this class strictly applies scoring
    logic without handling persistence or configuration loading.
    """
    
    def __init__(self, rules: RuleSet):
        """Initialize with a specific RuleSet."""
        self.rules = rules
    
    def _get_rule(self, name: str) -> Optional[RuleParameter]:
        return self.rules.rules.get(name)

    @staticmethod
    def _best_acceptor(core: str, params) -> Optional[Tuple[int, float]]:
        """Return best acceptor position (0-based) and weight using 1-based rule positions."""
        best_pos = None
        best_weight = None
        for pos in params.allowed_positions:
            idx = pos - 1
            if idx < len(core) and core[idx] in params.residues:
                weight = params.position_weights.get(str(pos), 1.0)
                if best_weight is None or weight > best_weight:
                    best_weight = weight
                    best_pos = idx
        if best_pos is None:
            return None
        return best_pos, best_weight

    def score_cleavage_site(
        self,
        sequence: str,
        position: int,
        orientation: Orientation = "leader_core",
    ) -> Tuple[float, Dict[str, float]]:
        """
        Score a potential cleavage site.
        
        Args:
            sequence: Full precursor sequence
            position: Cleavage position (0-indexed, core starts at this position)
        
        Returns:
            Tuple of (total_score, score_breakdown)
        """
        if orientation == "leader_core":
            leader = sequence[:position]
            core = sequence[position:]
        else:
            core = sequence[:position]
            leader = sequence[position:]
        scores = {}
        
        # Rule: Core start residue
        rule = self._get_rule("core_start_residue")
        if rule and rule.enabled and core:
            aa = core[0]
            params = rule.parameters
            if aa in params.hard_forbidden:
                scores["core_start"] = rule.weight * params.penalty_forbidden
            elif aa in params.scores:
                scores["core_start"] = rule.weight * params.scores[aa]
            else:
                scores["core_start"] = rule.weight * params.penalty_other
        
        # Rule: Ring closure residue (D/E at position 7-9)
        rule = self._get_rule("ring_closure_residue")
        if rule and rule.enabled:
            params = rule.parameters
            if params.allowed_positions and params.residues:
                best = self._best_acceptor(core, params)
                if best is not None:
                    _, best_weight = best
                    scores["ring_closure"] = rule.weight * best_weight
                else:
                    scores["ring_closure"] = rule.weight * params.penalty_if_missing_in_allowed_window
        
        # Rule: Core length
        rule = self._get_rule("core_length")
        if rule and rule.enabled:
            params = rule.parameters
            if params.optimal_min <= len(core) <= params.optimal_max:
                scores["core_length"] = rule.weight * params.score_optimal
            elif params.min <= len(core) <= params.max:
                scores["core_length"] = rule.weight * params.score_acceptable
            else:
                scores["core_length"] = rule.weight * params.penalty_outside
        
        # Rule: Leader length
        rule = self._get_rule("leader_length")
        if rule and rule.enabled:
            params = rule.parameters
            if params.min <= len(leader) <= params.max:
                scores["leader_length"] = rule.weight
            else:
                scores["leader_length"] = -rule.weight
        
        # Rule: Cleavage motif (residue before + first residue of core)
        rule = self._get_rule("cleavage_motif")
        if rule and rule.enabled and core and leader:
            if orientation == "leader_core" and position > 0:
                motif = sequence[position - 1] + core[0]
            else:
                motif = core[-1] + leader[0]
            motifs = rule.parameters.motifs
            if motif in motifs:
                scores["cleavage_motif"] = rule.weight * motifs[motif] / 2.0

        # Rule: Leader C-terminal motifs (only meaningful for leader_core orientation)
        rule = self._get_rule("leader_motifs")
        if rule and rule.enabled and len(leader) >= 3 and orientation == "leader_core":
            leader_suffix = leader[-6:]
            for motif in rule.parameters.motifs:
                if motif in leader_suffix:
                    scores["leader_motif"] = rule.weight
                    break

        # Rule: Leader penultimate residue (only meaningful for leader_core orientation)
        rule = self._get_rule("leader_penultimate_residue")
        if rule and rule.enabled and len(leader) >= 2 and orientation == "leader_core":
            penultimate = leader[-2]
            params = rule.parameters
            if penultimate in params.preferred:
                scores["leader_penultimate"] = rule.weight * params.preferred[penultimate]
            elif penultimate in params.allowed_rare:
                scores["leader_penultimate"] = rule.weight * params.allowed_rare[penultimate]
            else:
                scores["leader_penultimate"] = rule.weight * params.penalty_other

        # Rule: Leader signature motifs (optional regex)
        rule = self._get_rule("leader_signature_motifs")
        if rule and rule.enabled and leader:
            params = rule.parameters
            best = None
            if params.allow_regex:
                import re
                for pattern, weight in params.patterns.items():
                    if re.search(pattern, leader):
                        best = weight if best is None else max(best, weight)
            else:
                for pattern, weight in params.patterns.items():
                    if pattern in leader:
                        best = weight if best is None else max(best, weight)
            if best is not None:
                scores["leader_signature"] = rule.weight * best
        
        # Rule: Glycine-rich core
        rule = self._get_rule("glycine_rich_core")
        if rule and rule.enabled and core:
            gly_ratio = core.count("G") / len(core)
            if gly_ratio >= rule.parameters.min_ratio:
                scores["glycine_rich"] = rule.weight
        
        # Rule: Hydrophobicity transition (Kyte-Doolittle)
        rule = self._get_rule("hydrophobicity_transition")
        if rule and rule.enabled:
            window = rule.parameters.window_size or min(len(leader), len(core))
            if orientation == "leader_core":
                leader_window = leader[-window:] if window else leader
                core_window = core[:window] if window else core
            else:
                leader_window = leader[:window] if window else leader
                core_window = core[-window:] if window else core
            if leader_window and core_window:
                leader_hydro = sum(_KYTE_DOOLITTLE.get(aa, 0.0) for aa in leader_window) / len(leader_window)
                core_hydro = sum(_KYTE_DOOLITTLE.get(aa, 0.0) for aa in core_window) / len(core_window)
                diff = core_hydro - leader_hydro
                if diff >= rule.parameters.min_difference:
                    scores["hydrophobicity"] = rule.weight * min(diff / 2, 1.5)
        
        # Rule: Charge transition (net charge at pH 7.0)
        rule = self._get_rule("charge_transition")
        if rule and rule.enabled and leader and core:
            leader_charge = ProteinAnalysis(leader).charge_at_pH(7.0)
            core_charge = ProteinAnalysis(core).charge_at_pH(7.0)
            diff = leader_charge - core_charge if rule.parameters.leader_more_charged else core_charge - leader_charge
            if diff >= rule.parameters.min_difference:
                scores["charge_transition"] = rule.weight * min(diff * 2, 1.0)
        
        # Rule: Tail lock residue
        rule = self._get_rule("tail_lock_residue")
        if rule and rule.enabled:
            region_start = rule.parameters.region_start
            if len(core) > region_start:
                tail = core[region_start:]
                bulky_residues = set(rule.parameters.residues)
                if sum(aa in bulky_residues for aa in tail) >= rule.parameters.require_at_least:
                    scores["tail_lock"] = rule.weight

        # Rule: Steric lock region relative to acceptor position
        rule = self._get_rule("steric_lock_region")
        if rule and rule.enabled and core:
            params = rule.parameters
            acceptor_pos = None
            ring_rule = self._get_rule("ring_closure_residue")
            if ring_rule and ring_rule.enabled:
                ring_params = ring_rule.parameters
                best = self._best_acceptor(core, ring_params)
                acceptor_pos = best[0] if best is not None else None
            if acceptor_pos is not None:
                start_offset = params.relative_to_acceptor.get("start_offset", 0)
                end_offset = params.relative_to_acceptor.get("end_offset", 0)
                start = acceptor_pos + start_offset
                end = acceptor_pos + end_offset
                if start < 0:
                    start = 0
                if end >= len(core):
                    end = len(core) - 1
                if end >= start:
                    region = core[start:end + 1]
                    bulky = set(params.bulky_residues)
                    if sum(aa in bulky for aa in region) >= params.require_at_least:
                        scores["steric_lock"] = rule.weight
        
        # Rule: Forbidden core patterns
        rule = self._get_rule("forbidden_core_patterns")
        if rule and rule.enabled and core:
            for pattern in rule.parameters.patterns:
                if pattern in core:
                    scores["forbidden_pattern"] = rule.weight * rule.parameters.penalty_per_match
                    break  # Only penalize once
        
        total = float(sum(scores.values()))
        return total, scores


class CorePredictor:
    """
    Predict core sequences from lasso peptide precursors.
    """
    
    def __init__(self, rule_engine: Optional[RuleEngine] = None):
        """Initialize with a RuleEngine instance. If None, loads active rules."""
        if rule_engine is None:
            rules = load_active_ruleset()
            self.rules = RuleEngine(rules)
        else:
            self.rules = rule_engine
    
    def predict(
        self,
        sequence: str,
        top_n: int = 3,
        min_leader: int = 10,
        min_core: Optional[int] = None,
        allow_inverted: bool = False,
    ) -> PrecursorPrediction:
        """
        Predict core sequence(s) from a precursor.
        Optionally uses related precursors for conservation analysis.
        """
        sequence = sequence.strip().upper()
        predictions = []
        
        if min_core is None:
            core_rule = self.rules._get_rule("core_length")
            if core_rule and core_rule.enabled:
                min_core = core_rule.parameters.min
            else:
                min_core = 10

        def _iter_positions(seq_len: int, orientation: Orientation):
            if orientation == "leader_core":
                start = min_leader
                end = seq_len - min_core
            else:
                start = min_core
                end = seq_len - min_leader
            if end < start:
                return range(0)
            return range(start, end + 1)

        orientations: List[Orientation] = ["leader_core"]
        if allow_inverted:
            orientations.append("core_leader")

        for orientation in orientations:
            for pos in _iter_positions(len(sequence), orientation):
                score, breakdown = self.rules.score_cleavage_site(sequence, pos, orientation=orientation)

                if orientation == "leader_core":
                    leader = sequence[:pos]
                    core = sequence[pos:]
                else:
                    core = sequence[:pos]
                    leader = sequence[pos:]
            
                # Determine confidence
                if score >= 8.0:
                    confidence = "high"
                elif score >= 5.0:
                    confidence = "medium"
                else:
                    confidence = "low"
            
                # Include which rules were responsible for the score
                reasons = [
                    f"{name}:{value:+.2f}"
                    for name, value in sorted(
                        breakdown.items(),
                        key=lambda item: abs(item[1]),
                        reverse=True,
                    )
                    if value
                ]
            
                result = PredictionResult(
                    cleavage_site=pos,
                    leader=leader,
                    core=core,
                    score=score,
                    score_breakdown=breakdown,
                    confidence=confidence,
                    reasons=reasons,
                    orientation=orientation,
                )
                predictions.append(result)
        
        # Sort by score and take top N
        predictions.sort(key=lambda p: p.score, reverse=True)
        top_predictions = predictions[:top_n]

        return PrecursorPrediction(
            sequence=sequence,
            length=len(sequence),
            predictions=top_predictions,
            best_prediction=top_predictions[0] if top_predictions else None
        )
    
    def predict_batch(
        self,
        sequences: List[str],
        top_n: int = 3
    ) -> List[PrecursorPrediction]:
        """Predict cores for multiple sequences."""
        return [self.predict(seq, top_n) for seq in sequences]

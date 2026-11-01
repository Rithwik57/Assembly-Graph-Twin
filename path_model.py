"""Optional PyTorch path scorer for reference-labeled bubble candidates.

Install the optional dependency from requirements-ml.txt. The default tool
does not import PyTorch unless a model command or --model-checkpoint is used.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

FEATURE_COUNT = 6
MODEL_VERSION = "graph-path-scorer-0.1"


def _torch():
    try:
        import torch
        from torch import nn
    except ImportError as error:
        raise RuntimeError("PyTorch is optional; install requirements-ml.txt from the project folder") from error
    return torch, nn


def _read_examples(paths: list[Path]) -> list[dict[str, Any]]:
    examples = []
    for path in paths:
        report = json.loads(path.read_text(encoding="utf-8"))
        for item in report.get("training_examples", []):
            if len(item.get("features", [])) != FEATURE_COUNT or item.get("label") not in (0, 1):
                raise ValueError(f"{path}: malformed training example")
            examples.append(item)
    if not examples or len({item["label"] for item in examples}) < 2:
        raise ValueError("training reports must contain labeled positive and negative paths")
    return examples


def _model():
    torch, nn = _torch()

    class PathScorer(nn.Module):
        def __init__(self):
            super().__init__()
            self.layers = nn.Sequential(nn.Linear(FEATURE_COUNT, 16), nn.ReLU(), nn.Linear(16, 1))

        def forward(self, values):
            return self.layers(values).squeeze(-1)

    return torch, PathScorer


def train_model(train_reports: list[Path], validation_reports: list[Path], output: Path, seed: int = 1) -> dict[str, Any]:
    torch, PathScorer = _model()
    torch.manual_seed(seed)
    train_groups = {
        json.loads(path.read_text(encoding="utf-8")).get("inputs", {}).get("dataset_id")
        for path in train_reports
    }
    validation_groups = {
        json.loads(path.read_text(encoding="utf-8")).get("inputs", {}).get("dataset_id")
        for path in validation_reports
    }
    if None in train_groups or None in validation_groups:
        raise ValueError("all benchmark reports must include dataset_id; regenerate reports with explicit genome/sample IDs")
    if train_groups & validation_groups:
        raise ValueError("training and validation dataset_id values must be disjoint")
    train = _read_examples(train_reports)
    validation = _read_examples(validation_reports)
    train_x = torch.tensor([item["features"] for item in train], dtype=torch.float32)
    train_y = torch.tensor([item["label"] for item in train], dtype=torch.float32)
    mean = train_x.mean(dim=0)
    scale = train_x.std(dim=0, unbiased=False).clamp(min=1e-6)
    train_x = (train_x - mean) / scale
    model = PathScorer()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
    loss_fn = torch.nn.BCEWithLogitsLoss()
    for _ in range(250):
        optimizer.zero_grad()
        loss = loss_fn(model(train_x), train_y)
        loss.backward()
        optimizer.step()

    val_x = torch.tensor([item["features"] for item in validation], dtype=torch.float32)
    val_y = torch.tensor([item["label"] for item in validation], dtype=torch.float32)
    with torch.no_grad():
        logits = model((val_x - mean) / scale).tolist()
    # Temperature scaling minimizes validation Brier score on a deterministic grid.
    temperatures = [0.25 + index * 0.05 for index in range(96)]
    def sigmoid(value: float) -> float:
        return 1.0 / (1.0 + math.exp(-max(-40.0, min(40.0, value))))
    temperature = min(
        temperatures,
        key=lambda item: sum((sigmoid(logit / item) - float(label)) ** 2 for logit, label in zip(logits, val_y.tolist())) / len(logits),
    )
    probabilities = [sigmoid(logit / temperature) for logit in logits]
    # Select the largest validation coverage whose selected predictions have >=95% precision.
    threshold = 1.01
    for candidate_threshold in sorted(set(probabilities)):
        selected = [probability for probability in probabilities if probability >= candidate_threshold]
        labels = [int(label) for probability, label in zip(probabilities, val_y.tolist()) if probability >= candidate_threshold]
        if selected and sum(labels) / len(labels) >= 0.95:
            threshold = candidate_threshold
            break
    state = {name: value.detach().cpu().tolist() for name, value in model.state_dict().items()}
    checkpoint = {
        "model_version": MODEL_VERSION,
        "feature_count": FEATURE_COUNT,
        "feature_mean": mean.tolist(),
        "feature_scale": scale.tolist(),
        "temperature": temperature,
        "auto_resolve_threshold": threshold,
        "validation_examples": len(validation),
        "validation_brier_score": sum((probability - float(label)) ** 2 for probability, label in zip(probabilities, val_y.tolist())) / len(probabilities),
        "train_reports": [str(path) for path in train_reports],
        "validation_reports": [str(path) for path in validation_reports],
        "train_dataset_ids": sorted(group for group in train_groups if group),
        "validation_dataset_ids": sorted(group for group in validation_groups if group),
        "train_dataset_ids": sorted(group for group in train_groups if group),
        "validation_dataset_ids": sorted(group for group in validation_groups if group),
        "weights": state,
    }
    output.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {key: value for key, value in checkpoint.items() if key != "weights"}


def make_scorer(checkpoint: dict[str, Any]):
    torch, PathScorer = _model()
    model = PathScorer()
    state = {name: torch.tensor(value, dtype=torch.float32) for name, value in checkpoint["weights"].items()}
    model.load_state_dict(state)
    model.eval()

    def score(features: list[float]) -> float:
        values = (torch.tensor(features, dtype=torch.float32) - torch.tensor(checkpoint["feature_mean"])) / torch.tensor(checkpoint["feature_scale"])
        with torch.no_grad():
            logit = float(model(values).item()) / float(checkpoint["temperature"])
        return 1.0 / (1.0 + math.exp(-max(-40.0, min(40.0, logit))))

    return score


def score_feature_rows(feature_rows: dict[str, list[float]], checkpoint: dict[str, Any]) -> dict[str, float]:
    """Score all paths in a bubble using the optional neural model."""
    score = make_scorer(checkpoint)
    return {branch: score(features) for branch, features in feature_rows.items()}


def resolve_snapshot(snapshot: dict[str, Any], checkpoint_path: Path) -> dict[str, Any]:
    """Return per-candidate NN proposals; never changes or removes graph paths."""
    from benchmark import path_features

    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    if checkpoint.get("feature_count") != FEATURE_COUNT:
        raise ValueError("checkpoint feature schema is incompatible")
    candidates = []
    score = make_scorer(checkpoint)
    for diagnostic in snapshot["diagnostics"]:
        if diagnostic["type"] != "bubble":
            continue
        bubble = diagnostic["evidence"]
        probabilities = {branch: score(path_features(bubble, branch, snapshot["graph"])) for branch in bubble["branches"]}
        ranked = sorted(probabilities, key=lambda branch: (-probabilities[branch], branch))
        selected = ranked[0] if ranked and probabilities[ranked[0]] >= checkpoint["auto_resolve_threshold"] else None
        candidates.append({
            "candidate_id": bubble["candidate_id"],
            "selected_branch": selected,
            "decision": "auto_resolved" if selected is not None else "abstain",
            "branch_probabilities": probabilities,
            "threshold": checkpoint["auto_resolve_threshold"],
            "model_version": checkpoint["model_version"],
            "rationale": "highest calibrated path score exceeded the validation-selected threshold" if selected else "no path exceeded the validation-selected threshold",
            "alternatives_preserved": list(bubble["branches"]),
        })
    return {
        "model_version": checkpoint["model_version"],
        "checkpoint_sha256": hashlib.sha256(checkpoint_path.read_bytes()).hexdigest(),
        "candidates": candidates,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    train = sub.add_parser("train")
    train.add_argument("--train-reports", nargs="+", type=Path, required=True)
    train.add_argument("--validation-reports", nargs="+", type=Path, required=True)
    train.add_argument("--seed", type=int, default=1)
    train.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        print(json.dumps(train_model(args.train_reports, args.validation_reports, args.output, args.seed), indent=2))
    except (OSError, ValueError, RuntimeError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()

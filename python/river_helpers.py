"""Helpers for calling River regressors from R via reticulate."""

from __future__ import annotations

import importlib
import math
import numbers
from typing import Any, Iterable


def _coerce_scalar(value: Any) -> Any:
    if hasattr(value, "item"):
        try:
            value = value.item()
        except Exception:
            pass

    if isinstance(value, bool):
        return value
    if isinstance(value, numbers.Integral):
        return int(value)
    if isinstance(value, numbers.Real):
        return float(value)
    return value


def _coerce_feature_dict(row: Any) -> dict[str, Any]:
    if isinstance(row, dict):
        items = row.items()
    else:
        try:
            items = row.items()
        except AttributeError as exc:
            raise TypeError("Expected a mapping-like object for a River feature row.") from exc

    return {str(key): _coerce_scalar(value) for key, value in items}


def _resolve_attr(module_name: str, attr_name: str) -> Any:
    if not module_name.startswith("river."):
        module_name = f"river.{module_name}"
    module = importlib.import_module(module_name)
    return getattr(module, attr_name)


def _materialize(value: Any) -> Any:
    if isinstance(value, dict):
        if "class" in value and ("module" in value or "qualified_name" in value):
            if "qualified_name" in value:
                qualified_name = value["qualified_name"]
                module_name, attr_name = qualified_name.rsplit(".", 1)
            else:
                module_name = value["module"]
                attr_name = value["class"]

            params = value.get("params", {}) or {}
            params = {key: _materialize(val) for key, val in params.items()}
            cls = _resolve_attr(module_name, attr_name)
            return cls(**params)

        return {key: _materialize(val) for key, val in value.items()}

    if isinstance(value, (list, tuple)):
        return [_materialize(item) for item in value]

    return _coerce_scalar(value)


def _candidate_model_paths(model_name: str) -> list[str]:
    aliases = {
        "ARFRegressor": [
            "river.forest.ARFRegressor",
            "river.ensemble.AdaptiveRandomForestRegressor",
        ],
        "AdaptiveRandomForestRegressor": [
            "river.forest.ARFRegressor",
            "river.ensemble.AdaptiveRandomForestRegressor",
        ],
    }

    paths: list[str] = []
    if model_name in aliases:
        paths.extend(aliases[model_name])

    if "." in model_name:
        if model_name.startswith("river."):
            paths.append(model_name)
        else:
            paths.append(f"river.{model_name}")
    else:
        for module_name in (
            "forest",
            "tree",
            "ensemble",
            "linear_model",
            "neighbors",
            "dummy",
        ):
            paths.append(f"river.{module_name}.{model_name}")

    seen: set[str] = set()
    ordered_paths: list[str] = []
    for path in paths:
        if path not in seen:
            seen.add(path)
            ordered_paths.append(path)
    return ordered_paths


def create_model(model_name: str = "ARFRegressor", control: dict[str, Any] | None = None) -> Any:
    kwargs = _materialize(control or {})
    last_error: Exception | None = None

    for path in _candidate_model_paths(model_name):
        module_name, attr_name = path.rsplit(".", 1)
        try:
            cls = _resolve_attr(module_name, attr_name)
            return cls(**kwargs)
        except (ImportError, AttributeError, TypeError, ValueError) as exc:
            last_error = exc

    raise ValueError(
        f"Unable to resolve River regressor '{model_name}'. Last error: {last_error}"
    )


def fit_rows(model: Any, rows: Iterable[Any], y: Iterable[Any]) -> Any:
    for row, target in zip(rows, y):
        updated = model.learn_one(_coerce_feature_dict(row), _coerce_scalar(target))
        if updated is not None:
            model = updated
    return model


def predict_one(model: Any, row: Any) -> float:
    prediction = model.predict_one(_coerce_feature_dict(row))
    if prediction is None:
        return math.nan
    return float(_coerce_scalar(prediction))


def learn_one(model: Any, row: Any, y: Any) -> Any:
    updated = model.learn_one(_coerce_feature_dict(row), _coerce_scalar(y))
    return model if updated is None else updated

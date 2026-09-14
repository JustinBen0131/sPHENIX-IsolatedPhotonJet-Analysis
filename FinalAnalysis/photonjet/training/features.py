"""Ordered H70 photon-identification feature contracts."""

PP_H70_FEATURES = (
    "cluster_Et",
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "vertexz",
    "cluster_Eta",
    "e11_over_e33",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "e32_over_e35",
)

AUAU_H70_FEATURES = (
    "cluster_Et",
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_weta33_cogx",
    "cluster_wphi33_cogx",
    "vertexz",
    "cluster_Eta",
    "e11_over_e33",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "e32_over_e35",
    "centrality",
)


def feature_contract(system: str) -> tuple[str, ...]:
    if system == "pp":
        return PP_H70_FEATURES
    if system == "auau":
        return AUAU_H70_FEATURES
    raise ValueError("system must be pp or auau")

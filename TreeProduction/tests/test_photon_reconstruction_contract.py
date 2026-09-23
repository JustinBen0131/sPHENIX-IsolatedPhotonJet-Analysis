"""Bounded photon-source guards and extracted arithmetic checks.

These tests do not import or mock sPHENIX. One test compiles the actual standalone
arithmetic blocks from PhotonClusterBuilder.cc against ordinary arrays/scalars.
It does not compile the reconstruction module or certify ROOT/Fun4All behavior.
"""
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
BUILDER = (ROOT / "coresoftware/offline/packages/CaloReco/PhotonClusterBuilder.cc").read_text()
HEADER = (ROOT / "coresoftware/offline/packages/CaloReco/PhotonClusterBuilder.h").read_text()
CAPTURE = (ROOT / "TreeProduction/src/internal/Photons.cc").read_text()
PRODUCTION = (ROOT / "TreeProduction/src/internal/ProductionReconstruction.cc").read_text()
CONFIG = (ROOT / "TreeProduction/config/tree_production.yaml").read_text()


def body(source, signature):
    start = source.index("{", source.index(signature))
    depth = 1
    end = start + 1
    while depth:
        depth += (source[end] == "{") - (source[end] == "}")
        end += 1
    return source[start + 1:end - 1]


SHAPES = body(BUILDER, "bool PhotonClusterBuilder::calculate_shower_shapes(")
ISOLATION = body(BUILDER, "void PhotonClusterBuilder::calculate_isolation(")
EVENT = body(BUILDER, "int PhotonClusterBuilder::process_event(")
REGISTRATION = body(PRODUCTION, "void registerPhotons(")


class PhotonSourceContract(unittest.TestCase):
    def test_one_definition_per_builder_method(self):
        names = re.findall(r"^(?:int|void|bool|float|double|std::vector<int>) "
                           r"PhotonClusterBuilder::(\w+)\(", BUILDER, re.M)
        self.assertEqual(len(names), len(set(names)))
        for name in names:
            self.assertRegex(HEADER, rf"\b{name}\(")

    def test_same_configurable_h70_default_in_both_systems(self):
        for system in ("pp", "auau"):
            section = CONFIG.split(f"\n{system}:\n", 1)[1]
            section = re.split(r"\n[a-z_]+:\n", section, maxsplit=1)[0]
            self.assertRegex(section, r"shower_shape_tower_min_energy_gev:\s+0.070")
        self.assertIn("reco.showerShapeTowerMinEnergyGeV", REGISTRATION)
        self.assertIn("const double primaryFloor = namedShape", CAPTURE)
        self.assertIn('primary ? "" : "shower_" + definitionName + "_"', CAPTURE)
        self.assertIn("definitions: [H70, H0]", CONFIG)

    def test_explicit_independent_isolation_policy(self):
        self.assertIn("set_isolation_min_tower_energy(reco.isPP ? 0.12F : 0.0F)",
                      REGISTRATION)
        layer = body(BUILDER, "float PhotonClusterBuilder::calculate_layer_et(")
        self.assertIn("m_iso_min_tower_E > 0.0F && energy <= m_iso_min_tower_E", layer)
        self.assertIn("set_isolation_axis_cog(true, nativeAxisFloor)", REGISTRATION)
        self.assertIn("m_isolation_axis_floor", ISOLATION)

    def test_view_identity_uses_exact_stored_floor(self):
        self.assertIn("floor != static_cast<float>(reco.showerShapeTowerMinEnergyGeV)",
                      REGISTRATION)
        self.assertIn("primaryFloor == static_cast<double>(static_cast<float>(floor))",
                      CAPTURE)
        self.assertNotIn("std::abs(primaryFloor - floor)", CAPTURE)

    def test_h0_cannot_write_primary_ancillary_outputs(self):
        self.assertIn('"shower_" + view.first + "_", false', EVENT)
        self.assertNotIn('set_shower_shape_parameter("iso_', SHAPES)
        self.assertNotIn('set_shower_shape_parameter("ihcal_', SHAPES)
        timing = body(SHAPES, "if (include_ancillary)\n      ")
        self.assertIn("clusteravgtime +=", timing)
        self.assertIn('if (include_ancillary) photon->set_shower_shape_parameter("mean_time"', SHAPES)
        self.assertEqual(EVENT.count("calculate_isolation("), 1)
        self.assertEqual(EVENT.count("calculate_hcal_shapes("), 1)
        self.assertLess(SHAPES.index('"mean_time"'), SHAPES.index("showershape.size() < 12"))

    def test_nominal_kinematics_not_recomputed_in_capture(self):
        capture = body(CAPTURE, "void PhotonJetTree::capturePhotons(")
        self.assertNotIn("RawClusterUtility::", capture)
        for key in ("vertex_z", "cluster_pt", "cluster_eta", "cluster_phi"):
            self.assertIn(f'namedShape(cluster, "{key}")', capture)
            self.assertIn(f'set_shower_shape_parameter("{key}"', EVENT)
        self.assertIn("m_photonBuilder->get_vertex_z()", capture)
        witnesses = body(CAPTURE, "void PhotonJetTree::capturePhotonIsolationConstituents(")
        self.assertNotIn("m_event.recoVertexZ", witnesses)
        self.assertIn("photon.producerVertexZ", witnesses)

    def test_canonical_vertex_and_archived_boundary(self):
        self.assertIn("reco.isArchivedDoubleInteraction", REGISTRATION)
        self.assertIn("VertexSource::Mbd", REGISTRATION)
        self.assertIn("reco.isPP && reco.isSimulation", REGISTRATION)
        self.assertIn("VertexSource::GlobalMbd", REGISTRATION)
        self.assertIn("VertexSource::MbdThenGlobal", REGISTRATION)
        vertex = body(BUILDER, "float PhotonClusterBuilder::select_vertex_z(")
        self.assertIn("vertices->begin()->second", vertex)
        self.assertIn("find_vertexes(GlobalVertex::MBD)", vertex)
        self.assertNotIn("Truth", vertex)

    def test_topo_signed_double_skip_and_single_subtraction(self):
        topo = body(BUILDER, "bool PhotonClusterBuilder::calculate_topocluster_iso(")
        code = re.sub(r"//[^\n]*", "", topo)
        self.assertNotIn("->isValid(", code)
        self.assertIn("if (!topo) continue;", code)
        self.assertIn("if (!std::isfinite(topo_et)) continue;", code)
        self.assertIn("double sum = 0.0", code)
        self.assertIn("const double topo_eta", code)
        self.assertIn("const double topo_et", code)
        self.assertEqual(code.count("sum - candidate_et"), 1)
        self.assertNotIn("topo_et > 0", code)
        self.assertIn("result >= 1.0e8F", code)
        self.assertIn('"iso_topo_03_valid"', ISOLATION)
        self.assertIn('"iso_topo_04_valid"', ISOLATION)
        self.assertIn('"iso_topo_" + key + "_valid"', CAPTURE)

    def test_topo_reconstructed_before_photons(self):
        self.assertLess(REGISTRATION.index("server->registerSubsystem(topo)"),
                        REGISTRATION.index("server->registerSubsystem(photons)"))
        for call in ("set_noise(0.0053, 0.0351, 0.0684)",
                     "set_significance(4.0, 2.0, 1.0)",
                     "set_minE_local_max(1.0, 2.0, 0.5)",
                     "set_R_shower(0.025)", "set_absE(true)",
                     "set_use_only_good_towers(true)"):
            self.assertIn(call, REGISTRATION)

    def test_topo_default_axis_is_independent_of_layer_axis_switch(self):
        self.assertIn("m_isolation_axis_cog || m_do_topocluster_isolation", ISOLATION)
        self.assertIn("m_isolation_axis_cog ? cog_eta : cluster_eta", ISOLATION)
        self.assertIn('set_shower_shape_parameter("iso_topo_axis_eta", cog_eta)', ISOLATION)
        self.assertIn("calculate_topocluster_iso(cog_eta, cog_phi, ET", ISOLATION)

    def test_explicit_addressing_and_sub1_geometry(self):
        closest = body(BUILDER, "std::vector<int> PhotonClusterBuilder::find_closest_hcal_tower(")
        self.assertIn("TowerInfoDefs::encode_hcal(channel)", closest)
        self.assertIn("m_use_explicit_tower_channels", closest)
        layer = body(BUILDER, "float PhotonClusterBuilder::calculate_layer_et(")
        self.assertIn("TowerInfoDefs::encode_emcal(channel)", layer)
        self.assertIn("TowerInfoDefs::encode_hcal(channel)", layer)
        self.assertIn("m_emc_sub1_tower_container, m_geomIH, "
                      "RawTowerDefs::CalorimeterId::HCALIN", ISOLATION)
        self.assertNotIn("get_tower_at_key", CAPTURE)

    def test_models_do_not_select_candidates(self):
        self.assertIn("set_do_bdt(false)", REGISTRATION)
        self.assertNotIn("add_bdt_model(", REGISTRATION)
        model_tail = EVENT[EVENT.index("if (m_do_bdt)"):]
        self.assertNotIn("continue;", model_tail)
        self.assertIn("AddCluster(photon)", model_tail)
        self.assertIn("set_require_complete_shower_window(true, nativeAxisFloor)", REGISTRATION)

    def test_timing_population_and_sentinel(self):
        timing = body(SHAPES, "if (include_ancillary)\n      ")
        self.assertIn("get_time() * towerinfo->get_energy()", timing)
        self.assertIn("++cluster_time_tower_count", timing)
        self.assertNotIn("get_isGood", timing)
        self.assertNotIn("tower_floor", timing)
        self.assertIn("clusteravgtime = -999.0F", SHAPES)

    def test_all_captured_named_quantities_have_writers(self):
        writers = set(re.findall(r'(?:set_shape|set_shower_shape_parameter)\("([^"]+)"',
                                 BUILDER))
        reads = set(re.findall(r'value\("([^"]+)"\)', CAPTURE))
        reads |= set(re.findall(r'namedShape\(\s*cluster,\s*"([^"]+)"\s*\)', CAPTURE))
        self.assertEqual(reads - writers, set())
        for radius in ("03", "04"):
            for prefix in ("iso_", "iso_sub_"):
                for layer in ("emcal", "hcalin", "hcalout", "valid"):
                    self.assertIn(f"{prefix}{radius}_{layer}", writers)
            for suffix in ("", "_valid"):
                self.assertIn(f"iso_topo_{radius}{suffix}", writers)
        self.assertIn('value("moment33_valid") > 0.5', CAPTURE)

    @unittest.skipUnless(shutil.which("clang++"), "standard C++ compiler unavailable")
    def test_extracted_real_arithmetic(self):
        # Compile only actual standard-C++ arithmetic, not fake detector APIs.
        moments = SHAPES[SHAPES.index("  float e11 ="):SHAPES.index('  set_shape("et1"')]
        layer = body(BUILDER, "float PhotonClusterBuilder::calculate_layer_et(")
        floor = layer[layer.index("    if (!std::isfinite(energy))"):
                      layer.index("    float et =")]
        program = r'''
#include <array>
#include <cmath>
#include <limits>
#include <vector>
#include <cassert>
std::array<float, 8> moments(float E77[7][7], int E77_ownership[7][7],
                            float avg_eta=20.5F, float avg_phi=20.5F) {
  int maxieta=20, maxiphi=20;
''' + moments + r'''
  return {e11,e33,e32,e35,weta_cogx,wphi_cogx,weta33_cogx,wphi33_cogx};
}
std::vector<float> accepted(float m_iso_min_tower_E, float m_shape_min_tower_E,
                           bool subtracted, const std::vector<float>& energies) {
  std::vector<float> output;
  for (float energy : energies) {
''' + floor + r'''
    output.push_back(energy);
  }
  return output;
}
int main() {
  float e[7][7] = {}; int owned[7][7] = {};
  e[3][3]=10; owned[3][3]=1;
  e[3][4]=9;  // full-grid rectangle, never an owned-moment contribution
  e[4][3]=2; owned[4][3]=1;
  e[5][3]=4; owned[5][3]=1; // outside 3x3 denominator
  auto m=moments(e,owned);
  assert(m[0]==10 && m[1]==21);
  assert(m[4]==18.0F/16.0F && m[5]==0.0F);
  assert(m[6]==2.0F/12.0F && m[7]==0.0F);
  float center[7][7]={}; int one[7][7]={};
  center[3][3]=10; one[3][3]=1;
  m=moments(center,one,20.75F,20.75F);
  assert(m[4]==0 && m[5]==0 && m[6]==0 && m[7]==0);
  const float at=.12F, above=std::nextafter(at, 1.0F);
  const std::vector<float> energies={-2.0F,0.0F,.07F,at,above,
                                     std::numeric_limits<float>::quiet_NaN()};
  auto au=accepted(0.0F,.07F,true,energies);
  assert(au.size()==5 && au[0]==-2.0F && au[1]==0.0F);
  assert(au==accepted(0.0F,.5F,true,energies)); // no shower coupling
  assert(au==accepted(0.0F,.07F,false,energies)); // raw zero-floor also signed
  auto pp=accepted(.12F,0.0F,false,energies);
  assert(pp.size()==1 && pp[0]==above);
  assert(pp==accepted(.12F,.7F,false,energies));
}
'''
        with tempfile.TemporaryDirectory(prefix="photon-arithmetic-") as tmp:
            source = Path(tmp) / "arithmetic.cc"
            binary = Path(tmp) / "arithmetic"
            source.write_text(program)
            subprocess.run(["clang++", "-std=c++17", str(source), "-o", str(binary)],
                           check=True, capture_output=True, text=True)
            subprocess.run([str(binary)], check=True, capture_output=True, text=True)


if __name__ == "__main__":
    unittest.main()

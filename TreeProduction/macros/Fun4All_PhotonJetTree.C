/**
 * @file Fun4All_PhotonJetTree.C
 * @brief Native steering entry point for isolated-photon + jet eventTree production.
 *
 * This macro is the thin execution bridge between the Python production
 * orchestrator and the compiled TreeProduction implementation.
 *
 * Normal production enters through:
 *
 *   produce_trees.sh
 *       -> scripts/produce.py
 *           -> Fun4All_PhotonJetTree.C
 *               -> Production.cc / Setup.cc / PhotonJetTree
 *
 * By the time this macro is invoked, produce.py has already selected one exact
 * maintained source interval and supplied its source coordinates, input bundle,
 * configuration identity and relevant source-specific provenance.
 *
 * Responsibilities of this macro
 * ------------------------------
 *   1. Translate the launch arguments into a native production::Job.
 *   2. Resolve that Job plus the configured reconstruction profile into a fully
 *      validated production::Plan through LoadPlan().
 *   3. Expose the native output-mapping readiness report used by produce.py's
 *      preflight.
 *   4. For real execution, configure conditions, construct Fun4All, register
 *      reconstruction and input streams, and run the requested finite interval.
 *   5. Propagate event-loop or finalization failures so incomplete products
 *      cannot be mistaken for accepted TreeProduction output.
 *
 * Authority boundaries
 * --------------------
 * This file contains steering, not detector or analysis physics.
 *
 * It does not implement:
 *   - photon reconstruction,
 *   - jet reconstruction,
 *   - centrality,
 *   - isolation,
 *   - photon identification,
 *   - truth matching,
 *   - weighting,
 *   - event classification, or
 *   - eventTree serialization.
 *
 * Those behaviors live in the compiled TreeProduction library. Production.cc
 * resolves the typed native plan, Setup.cc binds conditions/reconstruction and
 * physical inputs, and PhotonJetTree plus its internal implementation owns the
 * event-level capture and output.
 *
 * Normal users should invoke produce_trees.sh rather than calling this macro
 * directly. Direct invocation is a low-level interface and requires the caller
 * to provide a self-consistent source identity and finite execution contract.
 */

#include <TreeProduction/PhotonJetTree.h>
#include <TreeProduction/internal/Production.h>

#include <fun4all/Fun4AllServer.h>

#include <TSystem.h>
#include <nlohmann/json.hpp>

#include <cstdint>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

R__LOAD_LIBRARY(libPhotonJetTree.so)

/**
 * Execute or preflight one resolved TreeProduction reconstruction piece.
 *
 * The argument list is the native launch contract serialized by produce.py.
 * Together, these values describe one exact source interval and bind that
 * interval to its population, configuration, input and optional recovery state.
 *
 * Two invocation modes exist:
 *
 *   checkOnly = true
 *     Resolve and validate the native Plan and report output-mapping readiness.
 *     No Fun4All server is created and no event is read.
 *
 *   checkOnly = false
 *     Execute the exact finite interval requested by produce.py.
 *
 * Parameters
 * ----------
 * profile
 *   Named reconstruction profile resolved by Production.cc.
 *
 * inputList
 *   Exact ordered physical-input bundle for this source piece.
 *
 * outputFile
 *   Candidate eventTree ROOT output for this invocation.
 *
 * run, segment
 *   Source run/segment coordinates carried into provenance and reconstruction.
 *
 * nEvents, firstEntry
 *   Exact finite source interval to process. firstEntry identifies the beginning
 *   of the interval in source order; nEvents is its size.
 *
 * configFile
 *   Original TreeProduction configuration whose exact bytes are provenance-bound.
 *
 * sourceFileOrdinal, sourceGlobalEntryBegin
 *   Source-order coordinates supplied by the maintained catalogue.
 *
 * inputFileSha256, sourceManifestSha256
 *   Source and catalogue provenance supplied by the production frontend.
 *
 * sample, period, siDiRole
 *   Production-population descriptors required by the native plan.
 *
 * verbosity
 *   Fun4All/producer diagnostic verbosity.
 *
 * checkOnly
 *   Request native plan/mapping validation without constructing the event loop.
 *
 * resolvedConfigurationFile
 *   Parsed JSON configuration envelope created by produce.py and bound to the
 *   original configuration digest.
 *
 * historicalRng
 *   Archived p+p double-interaction RNG replay state, empty when inapplicable.
 *
 * cemcRecoveryMap, cemcRecoveryMapSha256
 *   Run-bound Au+Au DATA CEMC status recovery payload, empty when inapplicable.
 *
 * Return value
 * ------------
 *   0 on successful preflight or completed production.
 *   1 on a native setup, event-loop, finalization or contract failure.
 */
int Fun4All_PhotonJetTree(const std::string& profile,
                          const std::string& inputList,
                          const std::string& outputFile,
                          const int run,
                          const int segment,
                          const std::int64_t nEvents = -1,
                          const std::int64_t firstEntry = 0,
                          const std::string& configFile = "config/tree_production.yaml",
                          const std::int64_t sourceFileOrdinal = -1,
                          const std::string& inputFileSha256 = "",
                          const std::string& sourceManifestSha256 = "",
                          const std::string& sample = "",
                          const std::string& period = "",
                          const std::string& siDiRole = "",
                          const int verbosity = 0,
                          const bool checkOnly = false,
                          const std::string& resolvedConfigurationFile = "",
                          const std::int64_t sourceGlobalEntryBegin = -1,
                          const std::string& historicalRng = "",
                          const std::string& cemcRecoveryMap = "",
                          const std::string& cemcRecoveryMapSha256 = "")
{
  using photonjet::production::Job;
  using photonjet::production::Plan;

  try
  {
    /*
     * ------------------------------------------------------------------------
     * Construct the native per-piece Job
     * ------------------------------------------------------------------------
     *
     * The Python driver has already selected the logical source and interval.
     * This layer performs no population discovery. It transfers the resolved
     * launch contract into the typed structure consumed by Production.cc.
     */
    Job job;
    job.profile = profile;
    job.inputList = inputList;
    job.outputFile = outputFile;
    job.run = run;
    job.segment = segment;
    job.firstEntry = firstEntry;
    job.nEvents = nEvents;
    job.sourceFileOrdinal = sourceFileOrdinal;
    job.sourceGlobalEntryBegin = sourceGlobalEntryBegin;

    /*
     * Bind the exact ordered input-list bytes used by this native invocation.
     *
     * `inputUriSha256` proves the actual list presented to the macro, while the
     * catalogue-supplied source/file/manifest digests preserve the higher-level
     * scientific identity of that source.
     *
     * This is especially important for schema-v2 catalogues, where produce.py
     * may materialize an embedded input bundle at a temporary filesystem path:
     * the temporary locator can change, but the exact bytes remain provable.
     */
    // Bind the actual ordered DST locations, alongside the catalogue's
    // content digest. Never leave the native source URI identity empty.
    job.inputUriSha256 = photonjet::production::detail::fileSha256(inputList);
    job.inputFileSha256 = inputFileSha256;
    job.sourceManifestSha256 = sourceManifestSha256;
    job.sample = sample;
    job.period = period;
    job.siDiRole = siDiRole;
    job.historicalRng = historicalRng;
    job.cemcRecoveryMap = cemcRecoveryMap;
    job.cemcRecoveryMapSha256 = cemcRecoveryMapSha256;
    job.verbosity = verbosity;
    job.resolvedConfigurationFile = resolvedConfigurationFile;

    /*
     * Hash the steering macro that is actually executing.
     *
     * The digest becomes part of production provenance and allows downstream
     * validation to distinguish outputs made by different steering revisions
     * even when all nominal runtime arguments are otherwise identical.
     */
    job.macroSha256 = photonjet::production::detail::fileSha256(__FILE__);

    /*
     * ------------------------------------------------------------------------
     * Resolve the complete native production Plan
     * ------------------------------------------------------------------------
     *
     * LoadPlan() is the authority boundary between orchestration and physics.
     * It combines the source-specific Job with the requested configuration
     * profile and validates the complete typed production state.
     *
     * The macro does not reinterpret profile values after this point.
     */
    Plan plan = photonjet::production::LoadPlan(configFile, job);

    /*
     * ------------------------------------------------------------------------
     * Native output-mapping preflight
     * ------------------------------------------------------------------------
     *
     * Query the actually loaded producer implementation for schema coverage.
     * produce.py searches stdout for the TREEPRODUCTION_MAPPING_JSON prefix and
     * verifies both readiness and schema identity before admitting production.
     *
     * A missing native filling rule is therefore discovered before the event
     * loop rather than after an expensive reconstruction has produced an
     * incomplete tree.
     */
    const auto mapping = photonjet::production::OutputMappingStatusJson(plan.tree);
    std::cout << "TREEPRODUCTION_MAPPING_JSON=" << mapping << std::endl;
    if (!nlohmann::json::parse(mapping).at("ready").get<bool>())
      throw std::runtime_error("native output mapping is incomplete; see mapping report");

    /*
     * checkOnly is the read-only native preflight path used by produce.py.
     *
     * Reaching this point proves that the requested configuration can be
     * resolved by the current native implementation and that its eventTree
     * mapping is complete. No server, reconstruction modules, inputs or output
     * are created in this mode.
     */
    if (checkOnly) return 0;

    /*
     * Native production must always have a finite event budget.
     *
     * The Python orchestrator determines piece boundaries. Allowing native
     * execution with an unspecified budget would bypass that bounded-source
     * contract.
     */
    if (nEvents < 0) throw std::runtime_error("native execution requires an exact finite source event budget");

    /*
     * ------------------------------------------------------------------------
     * Empty logical source
     * ------------------------------------------------------------------------
     *
     * Fun4All interprets run(0) as unbounded processing rather than "zero
     * events." An explicitly empty maintained source therefore takes a separate
     * path: instantiate the producer directly, book the complete schema and
     * provenance, then finalize it without constructing Fun4All or reading DSTs.
     *
     * This preserves a valid zero-entry product without risking accidental
     * unbounded input processing.
     */
    if (nEvents == 0)
    {
      // Fun4All run(0) means unbounded processing. Book/end the producer
      // explicitly instead, without creating a server or reading a DST.
      PhotonJetTree empty(plan.tree, "PhotonJetTree");
      if (empty.Init(nullptr) != 0) return 1;
      return empty.End(nullptr) == 0 ? 0 : 1;
    }

    /*
     * ------------------------------------------------------------------------
     * Bind run-dependent conditions
     * ------------------------------------------------------------------------
     *
     * ConfigureConditions() resolves and validates the exact calibration state
     * required by the Plan before any reconstruction module is allowed to run.
     *
     * Print() exposes the resolved native plan for production diagnostics and
     * provenance inspection.
     */
    photonjet::production::ConfigureConditions(plan);
    photonjet::production::Print(plan);

    /*
     * ------------------------------------------------------------------------
     * Create the Fun4All execution environment
     * ------------------------------------------------------------------------
     */
    Fun4AllServer* server = Fun4AllServer::instance();
    server->Verbosity(verbosity);

    /*
     * PhotonJetTree is the canonical eventTree producer.
     *
     * It owns the event-level capture lifecycle and receives the already
     * resolved tree configuration. The macro does not duplicate any of its
     * event, object, truth, labeling, weighting or serialization logic.
     */
    auto* producer = new PhotonJetTree(plan.tree, "PhotonJetTree");

    /*
     * ------------------------------------------------------------------------
     * Register reconstruction first, then physical inputs
     * ------------------------------------------------------------------------
     *
     * RegisterReconstruction() instantiates the exact detector/reconstruction
     * chain required by the resolved Plan and binds the PhotonJetTree producer
     * to the reconstruction objects it must observe.
     *
     * RegisterInputs() then attaches the resolved source streams with the
     * synchronization policy determined by the Plan.
     */
    photonjet::production::RegisterReconstruction(server, plan, producer);
    photonjet::production::RegisterInputs(server, plan);
    photonjet::production::WriteProgress(plan.tree, "initialized", 0, 0, "", true);

    /*
     * ------------------------------------------------------------------------
     * Deterministic nonzero source offset
     * ------------------------------------------------------------------------
     *
     * Piece boundaries are source-entry boundaries, not independent datasets.
     * When this piece begins after entry zero, every physical input stream must
     * be advanced consistently to the requested firstEntry before the ordinary
     * event loop begins.
     *
     * StartAtFirstEntry() owns that synchronized positioning and verifies the
     * resulting source entry rather than relying on an unobserved skip.
     */
    if (firstEntry > 0)
    {
      // Open every source stream at the first entry; each event's entry is checked.
      photonjet::production::StartAtFirstEntry(server, plan);
      photonjet::production::WriteProgress(plan.tree, "deterministic_skip", 0, firstEntry, "", true);
    }

    /*
     * ------------------------------------------------------------------------
     * Execute the exact bounded reconstruction interval
     * ------------------------------------------------------------------------
     *
     * The event budget is supplied by produce.py's piece scheduler. Every event
     * reaching PhotonJetTree follows the single native capture spine defined by
     * PhotonJetTree::process_event().
     */
    const int status = server->run(static_cast<int>(nEvents));

    /*
     * A stopped/failed event loop must never acquire normal completion metadata.
     *
     * Mark the producer aborted before server->End() so finalization can close
     * the physical file while preserving the distinction between an interrupted
     * artifact and a certified completed product.
     */
    if (status != 0)
    {
      std::cerr << "Fun4All_PhotonJetTree: the event loop stopped with status " << status
                << "; the output is marked aborted" << std::endl;
      producer->markAborted();
    }

    /*
     * ------------------------------------------------------------------------
     * Finalize the Fun4All lifecycle
     * ------------------------------------------------------------------------
     *
     * PhotonJetTree::End() resolves terminal source accounting, writes completion
     * metadata and closes the output. The server is deleted only after every
     * registered subsystem has received End().
     */
    const int endStatus = server->End();
    delete server;

    /*
     * Either event-loop failure or End() failure invalidates the native
     * invocation. gSystem->Exit(1) ensures ROOT itself returns failure to
     * produce.py, which consequently refuses piece certification/publication.
     */
    if (status != 0 || endStatus != 0)
    {
      gSystem->Exit(1);
      return 1;
    }

    /*
     * Successful return means native execution and finalization completed.
     *
     * It does NOT by itself publish or accept the piece. produce.py subsequently
     * invokes libEventTree to verify completion metadata, source identity and
     * structural schema integrity before the file receives its checked receipt.
     */
    std::cout << "Fun4All_PhotonJetTree: complete" << std::endl;
    return 0;
  }
  catch (const std::exception& error)
  {
    /*
     * Any native contract, configuration, setup or reconstruction exception is
     * terminal for this invocation. Returning a ROOT process failure prevents
     * the Python orchestration layer from treating the candidate as complete.
     */
    std::cerr << "Fun4All_PhotonJetTree failed: " << error.what() << std::endl;
    gSystem->Exit(1);
    return 1;
  }
}

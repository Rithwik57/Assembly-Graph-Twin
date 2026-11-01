import unittest
import json
import importlib.util
import subprocess
import sys
from io import StringIO
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from assembly_graph_twin import (
    AssemblyGraphTwin,
    AssemblyGraph,
    Read,
    build_graph,
    classify_diagnostics,
    detect_bubbles,
    detect_disconnected_components,
    detect_tips,
    extract_kmers,
    graph_metrics,
    initial_reads,
    bridge_reads,
    parse_fasta,
    parse_fastq,
    reachable_within,
    reverse_complement,
)
from benchmark import _sample_rank, reference_supported_alignments, run_benchmark


def graph_from_edges(edges: list[tuple[str, str]]) -> AssemblyGraph:
    graph = AssemblyGraph(k=2)
    for index, (source, target) in enumerate(edges):
        graph.add_edge(source, target, source + target, f"edge-{index}")
    return graph


class AssemblyGraphTwinTests(unittest.TestCase):
    def test_extract_kmers(self):
        self.assertEqual(extract_kmers("ACTG", 3), ["ACT", "CTG"])

    def test_graph_tracks_edges_and_read_support(self):
        graph = build_graph([Read("r1", "ACTG"), Read("r2", "ACTG")], 3)
        self.assertEqual(set(graph.nodes), {"AC", "CT", "TG", "CA", "AG", "GT"})
        self.assertEqual(graph.edges[("AC", "CT")].supporting_reads, {"r1", "r2"})
        self.assertEqual(graph.edges[("AG", "GT")].supporting_reads, {"r1", "r2"})
        self.assertEqual(graph.nodes["AC"].coverage, 2)

    def test_bidirected_graph_is_invariant_to_read_orientation(self):
        sequence = "AACGTTGCA"
        reverse = reverse_complement(sequence)
        forward_graph = build_graph([Read("forward", sequence)], 4)
        reverse_graph = build_graph([Read("reverse", reverse)], 4)
        edge_support = lambda graph: {
            key: edge.support_count for key, edge in graph.edges.items()
        }
        self.assertEqual(graph_metrics(forward_graph), graph_metrics(reverse_graph))
        self.assertEqual(edge_support(forward_graph), edge_support(reverse_graph))
        compact_forward = build_graph([Read("forward", sequence)], 4, compact=True)
        compact_reverse = build_graph([Read("reverse", reverse)], 4, compact=True)
        self.assertEqual(graph_metrics(compact_forward), graph_metrics(compact_reverse))
        self.assertEqual(
            {key: edge.support_count for key, edge in compact_forward.edges.items()},
            {key: edge.support_count for key, edge in compact_reverse.edges.items()},
        )

    def test_bidirected_edges_are_only_observed_edges_and_valid_reverse_complement_twins(self):
        sequence = "AACG"
        k = 3
        graph = build_graph([Read("r1", sequence)], k)
        expected = set()
        for oriented in (sequence, reverse_complement(sequence)):
            expected.update(
                (oriented[index:index + k - 1], oriented[index + 1:index + k])
                for index in range(len(oriented) - k + 1)
            )
        self.assertEqual(set(graph.edges), expected)
        self.assertNotIn(("AA", "CG"), graph.edges)
        self.assertEqual(graph.strand_mode, "bidirected")

    def test_bidirected_mode_connects_overlapping_reads_from_opposite_strands(self):
        genome_segment = "AAACGTCAGT"
        reads = [
            Read("forward", genome_segment[:7]),
            Read("reverse", reverse_complement(genome_segment[4:])),
        ]
        directed = build_graph(reads, 4, strand_mode="forward")
        bidirected = build_graph(reads, 4, strand_mode="bidirected")
        compact = build_graph(reads, 4, compact=True, strand_mode="bidirected")
        self.assertEqual(graph_metrics(directed)["component_count"], 2)
        self.assertEqual(graph_metrics(bidirected)["component_count"], 1)
        self.assertEqual(graph_metrics(compact), graph_metrics(bidirected))
        self.assertLess(
            graph_metrics(bidirected)["strand_collapsed_node_count"],
            graph_metrics(bidirected)["node_count"],
        )

    def test_compact_graph_preserves_counts_and_topology_without_read_id_sets(self):
        reads = [
            Read("ref-1", "AAACCTG"),
            Read("ref-2", "AAACCTG"),
            Read("alt-1", "AAATCTG"),
            Read("repeat", "CCCC"),
        ]
        full_graph = build_graph(reads, 3)
        compact_graph = build_graph(reads, 3, compact=True)
        self.assertEqual(graph_metrics(compact_graph), graph_metrics(full_graph))
        self.assertEqual(compact_graph.read_node_paths, {})
        self.assertTrue(all(node.supporting_reads is None for node in compact_graph.nodes.values()))
        self.assertTrue(all(edge.supporting_reads is None for edge in compact_graph.edges.values()))
        self.assertEqual(
            {edge_key: edge.support_count for edge_key, edge in compact_graph.edges.items()},
            {edge_key: len(edge.supporting_reads) for edge_key, edge in full_graph.edges.items()},
        )
        self.assertEqual(
            {node_key: node.coverage for node_key, node in compact_graph.nodes.items()},
            {node_key: node.coverage for node_key, node in full_graph.nodes.items()},
        )
        self.assertEqual(compact_graph.edges[("CC", "CC")].support_count, 1)
        self.assertAlmostEqual(compact_graph.edges[("CC", "CC")].quality_support, 1.0)
        self.assertEqual(len(detect_bubbles(compact_graph)), len(detect_bubbles(full_graph)))
        with self.assertRaisesRegex(ValueError, "unique read IDs"):
            build_graph([Read("same", "AAAA"), Read("same", "CCCC")], 3, compact=True)

    def test_demo_contains_tip_and_bubble(self):
        graph = build_graph(initial_reads(), 3)
        self.assertTrue(detect_tips(graph))
        self.assertTrue(detect_bubbles(graph))

    def test_reachable_within_searches_multiple_hops(self):
        graph = build_graph([Read("left", "ABDE"), Read("right", "ACFE")], 2, strand_mode="forward")
        self.assertEqual(reachable_within(graph, "B", depth=3), {"D", "E"})
        self.assertEqual(reachable_within(graph, "C", depth=3), {"F", "E"})
        bubbles = detect_bubbles(graph)
        self.assertEqual(len(bubbles), 1)
        self.assertEqual(bubbles[0]["rejoin"], "E")

    def test_bubble_reports_nearest_rejoin_not_lexically_first_later_node(self):
        graph = graph_from_edges([
            ("S", "B"), ("S", "C"), ("B", "Z"), ("Z", "A"),
            ("C", "Y"), ("Y", "Z"),
        ])
        bubbles = detect_bubbles(graph)
        self.assertEqual(len(bubbles), 1)
        self.assertEqual(bubbles[0]["rejoin"], "Z")

    def test_long_bubble_at_search_limit_is_detected(self):
        graph = graph_from_edges([
            ("S", "A"), ("S", "B"), ("A", "C"), ("C", "R"),
            ("B", "D"), ("D", "R"),
        ])
        bubbles = detect_bubbles(graph)
        self.assertEqual([(bubble["start"], bubble["rejoin"]) for bubble in bubbles], [("S", "R")])

    def test_bubble_beyond_search_limit_is_not_reported(self):
        graph = graph_from_edges([
            ("S", "A"), ("S", "B"), ("A", "C"), ("C", "D"),
            ("D", "E"), ("E", "R"), ("B", "F"), ("F", "G"),
            ("G", "H"), ("H", "R"),
        ])
        self.assertEqual(detect_bubbles(graph), [])

    def test_three_branches_are_grouped_into_one_bubble(self):
        graph = graph_from_edges([
            ("S", "A"), ("S", "B"), ("S", "C"),
            ("A", "R"), ("B", "R"), ("C", "R"),
        ])
        bubbles = detect_bubbles(graph)
        self.assertEqual(len(bubbles), 1)
        self.assertEqual(bubbles[0]["branches"], ["A", "B", "C"])
        self.assertEqual(set(bubbles[0]["paths"]), {"A", "B", "C"})
        for path in bubbles[0]["paths"].values():
            self.assertEqual(path[0], "S")
            self.assertEqual(path[-1], "R")
            self.assertTrue(
                all(
                    (path[index], path[index + 1]) in graph.edges
                    for index in range(len(path) - 1)
                )
            )

    def test_nested_bubbles_are_reported_at_their_distinct_divergences(self):
        graph = graph_from_edges([
            ("S", "A"), ("S", "B"), ("A", "X"), ("A", "Y"),
            ("X", "R"), ("Y", "R"), ("B", "R"),
        ])
        bubbles = detect_bubbles(graph)
        self.assertEqual(
            [(bubble["start"], bubble["rejoin"]) for bubble in bubbles],
            [("A", "R"), ("S", "R")],
        )

    def test_cycles_without_reconvergence_are_ignored(self):
        graph = graph_from_edges([
            ("S", "A"), ("S", "B"), ("A", "X"), ("X", "A"), ("B", "Y"),
        ])
        self.assertEqual(detect_bubbles(graph), [])

    def test_cycle_that_reaches_reconvergence_is_reported(self):
        graph = graph_from_edges([
            ("S", "A"), ("S", "B"), ("A", "X"), ("X", "A"),
            ("X", "R"), ("B", "Y"), ("Y", "R"),
        ])
        bubbles = detect_bubbles(graph)
        self.assertEqual(
            [(bubble["start"], bubble["rejoin"]) for bubble in bubbles],
            [("S", "R")],
        )

    def test_crossing_branches_do_not_create_later_false_rejoin(self):
        graph = graph_from_edges([
            ("S", "A"), ("S", "B"), ("A", "X"), ("X", "R"),
            ("B", "X"), ("B", "Y"), ("Y", "R"),
        ])
        bubbles = detect_bubbles(graph)
        self.assertEqual(
            [(bubble["start"], bubble["rejoin"]) for bubble in bubbles],
            [("B", "R"), ("S", "X")],
        )

    def test_graph_insertion_order_does_not_change_bubbles(self):
        edges = [
            ("S", "A"), ("S", "B"), ("A", "R"), ("B", "R"),
        ]
        first = detect_bubbles(graph_from_edges(edges))
        second = detect_bubbles(graph_from_edges(list(reversed(edges))))
        self.assertEqual(first, second)

    def test_candidate_ids_and_nested_relationships_are_stable(self):
        edges = [
            ("S", "A"), ("S", "B"), ("A", "X"), ("A", "Y"),
            ("X", "R"), ("Y", "R"), ("B", "R"),
        ]
        first = detect_bubbles(graph_from_edges(edges))
        second = detect_bubbles(graph_from_edges(list(reversed(edges))))
        self.assertEqual(
            [bubble["candidate_id"] for bubble in first],
            [bubble["candidate_id"] for bubble in second],
        )
        self.assertTrue(all(bubble["candidate_id"].startswith("bubble-") for bubble in first))
        outer, inner = next(bubble for bubble in first if bubble["start"] == "S"), next(
            bubble for bubble in first if bubble["start"] == "A"
        )
        outer_relation = outer["relationships"][0]
        inner_relation = inner["relationships"][0]
        self.assertEqual(outer_relation["candidate_id"], inner["candidate_id"])
        self.assertEqual(inner_relation["candidate_id"], outer["candidate_id"])
        self.assertEqual(outer_relation["type"], "nested")
        self.assertEqual(outer_relation["shared_internal_nodes"], ["X"])
        self.assertEqual(outer_relation["shared_edges"], ["A->X", "X->R"])
        self.assertEqual(
            len({relation["candidate_id"] for bubble in first for relation in bubble["relationships"]}),
            2,
        )

    def test_unrelated_candidates_have_no_relationships(self):
        graph = graph_from_edges([
            ("S", "A"), ("S", "B"), ("A", "R"), ("B", "R"),
            ("T", "C"), ("T", "D"), ("C", "U"), ("D", "U"),
        ])
        bubbles = detect_bubbles(graph)
        self.assertEqual(len(bubbles), 2)
        self.assertTrue(all(bubble["relationships"] == [] for bubble in bubbles))

    def test_configurable_depth_and_search_limit_are_explicit(self):
        graph = graph_from_edges([
            ("S", "A"), ("S", "B"), ("A", "C"), ("C", "D"),
            ("D", "E"), ("E", "R"), ("B", "F"), ("F", "G"),
            ("G", "H"), ("H", "R"),
        ])
        self.assertEqual(detect_bubbles(graph), [])
        self.assertEqual(len(detect_bubbles(graph, depth=4)), 1)
        limited = classify_diagnostics(
            graph, bubble_depth=4, max_bubble_search_nodes=3
        )
        limited_evidence = next(
            item["evidence"]
            for item in limited
            if item["type"] == "bubble_search_limited"
        )
        self.assertEqual(limited_evidence["max_nodes"], 3)
        self.assertTrue(any(limited_evidence["branch_search_limited"].values()))
        self.assertTrue(any(limited_evidence["branch_limit_reason"].values()))
        self.assertIn("absence of a bubble is not established", limited_evidence["reason"])

    def test_invalid_bubble_depth_and_cli_error(self):
        with self.assertRaisesRegex(ValueError, "between 1 and 12"):
            AssemblyGraphTwin(bubble_depth=0)
        with self.assertRaisesRegex(ValueError, "between 1 and 12"):
            classify_diagnostics(AssemblyGraph(k=2), bubble_depth=13)
        with self.assertRaisesRegex(ValueError, "between 2 and 256"):
            AssemblyGraphTwin(max_bubble_search_nodes=257)
        command = subprocess.run(
            [
                sys.executable,
                str(PROJECT_ROOT / "assembly_graph_twin.py"),
                "demo",
                "--bubble-depth",
                "0",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(command.returncode, 2)
        self.assertIn("bubble depth must be between 1 and 12", command.stderr)

    def test_parse_multiline_fasta(self):
        reads = parse_fasta(
            StringIO(">read-1 description\nAAAC\nGT\n>read-2\nTTAA\n"),
            Path("fixture.fa"),
        )
        self.assertEqual([(read.read_id, read.sequence) for read in reads], [
            ("read-1", "AAACGT"), ("read-2", "TTAA")
        ])

    def test_parse_fastq_converts_phred33_to_mean_accuracy(self):
        reads = parse_fastq(
            StringIO("@read-1 description\nACGT\n+\nIIII\n"),
            Path("fixture.fq"),
        )
        self.assertEqual(reads[0].read_id, "read-1")
        self.assertEqual(reads[0].sequence, "ACGT")
        self.assertAlmostEqual(reads[0].quality, 1.0 - 10 ** (-40 / 10))

    def test_bubble_includes_full_path_read_support(self):
        graph = build_graph([Read("left", "ABDE"), Read("right", "ACFE")], 2)
        bubble = detect_bubbles(graph)[0]
        self.assertEqual(bubble["path_read_support"]["B"], ["left"])
        self.assertEqual(bubble["path_read_support"]["C"], ["right"])

    def test_edge_only_support_is_not_full_path_support(self):
        graph = build_graph(
            [
                Read("edge-only", "SA"),
                Read("left-full", "SAXR"),
                Read("right-full", "SBYR"),
            ],
            2,
            strand_mode="forward",
        )
        bubble = detect_bubbles(graph)[0]
        self.assertEqual(bubble["branch_support"]["A"], 2)
        self.assertEqual(bubble["path_read_support"]["A"], ["left-full"])
        self.assertNotIn("edge-only", bubble["path_read_support"]["A"])

    def test_fastq_quality_is_recorded_on_supported_edges(self):
        reads = parse_fastq(
            StringIO("@q\nACGT\n+\nIIII\n"), Path("quality.fq")
        )
        graph = build_graph(reads, 2)
        self.assertAlmostEqual(graph.edges[("A", "C")].quality_support, reads[0].quality)

    def test_empty_graph_and_invalid_k(self):
        self.assertEqual(build_graph([], 3).nodes, {})
        self.assertEqual(build_graph([Read("short", "A")], 3).nodes, {})
        with self.assertRaises(ValueError):
            AssemblyGraphTwin(k=1)

    def test_bubble_detector_rejects_self_loop_false_positive(self):
        graph = build_graph([Read("loop", "AAA"), Read("branch", "AAB")], 2)
        self.assertEqual(detect_bubbles(graph), [])

    def test_disconnected_components_are_diagnostics(self):
        graph = build_graph([Read("first", "AAAA"), Read("second", "CCCC")], 3, strand_mode="forward")
        diagnostics = detect_disconnected_components(graph)
        self.assertEqual([item["nodes"] for item in diagnostics], [["AA"], ["CC"]])
        classified = classify_diagnostics(graph)
        self.assertEqual(
            sum(item["type"] == "disconnected_component" for item in classified),
            2,
        )

    def test_twin_records_two_snapshots(self):
        twin = AssemblyGraphTwin(k=3)
        first = twin.update(initial_reads(), "initial")
        first_before_second_update = json.dumps(first, sort_keys=True)
        second = twin.update(initial_reads() + bridge_reads(), "updated")
        self.assertEqual(len(twin.snapshots), 2)
        self.assertEqual(first["snapshot_id"], 1)
        self.assertEqual(second["snapshot_id"], 2)
        self.assertNotEqual(first["input_sha256"], second["input_sha256"])
        self.assertTrue(all("heuristic_score" in item for item in first["diagnostics"]))
        self.assertTrue(all("confidence" not in item for item in first["diagnostics"]))
        first_bubble = next(item for item in first["diagnostics"] if item["type"] == "bubble")
        second_bubble = next(item for item in second["diagnostics"] if item["type"] == "bubble")
        self.assertEqual(first_bubble["evidence"]["search_depth"], 3)
        self.assertIn("branch_support", second_bubble["evidence"])
        self.assertEqual(json.dumps(first, sort_keys=True), first_before_second_update)
        self.assertEqual(second["batch_read_count"], 2)
        self.assertTrue(second["state_delta"]["edge_support_increases"])

    def test_reference_benchmark_labels_unique_path_and_reports_abstention(self):
        reads = PROJECT_ROOT / "fixtures" / "known_bubble" / "reads.fasta"
        reference = PROJECT_ROOT / "fixtures" / "known_bubble" / "reference.fasta"
        report = run_benchmark([reads], reference, k=3, seed=7, max_reads=None)
        sample_one = run_benchmark([reads], reference, k=3, seed=11, max_reads=2)
        sample_two = run_benchmark([reads], reference, k=3, seed=11, max_reads=2)
        self.assertEqual(report["metrics"]["candidate_count"], 1)
        self.assertEqual(report["metrics"]["uniquely_labeled_candidate_count"], 1)
        self.assertEqual(report["metrics"]["accuracy_on_resolved_labeled"], 1.0)
        self.assertEqual(report["candidates"][0]["reference_label_status"], "unique_reference_path")
        self.assertTrue(report["training_examples"])
        self.assertEqual(sample_one["inputs"]["total_input_read_count"], 3)
        self.assertEqual(sample_one["inputs"]["selected_read_count"], 2)
        self.assertEqual(sample_one["inputs"]["read_fingerprint_sha256"], sample_two["inputs"]["read_fingerprint_sha256"])
        summary = report["inputs"]["diagnostic_summary"]
        self.assertEqual(summary["bubble_candidate_count"], report["metrics"]["candidate_count"])
        self.assertIn("largest_components", summary["examples"])
        self.assertEqual(report["protocol"]["protocol_version"], "3.0")
        self.assertEqual(report["protocol"]["strand_mode"], "bidirected")
        self.assertEqual(report["protocol"]["label_semantics"], "reference_compatibility_proxy")
        self.assertEqual(report["protocol"]["sampling_algorithm"], "sha256_bottom_k_by_read_id_v1")
        self.assertIn("bubble_search_limit_reasons", summary)

    def test_reference_labeler_aligns_a_path_with_a_substitution_and_reverse_complement(self):
        query = "ACGTTGCATCGATGCATGACTGACCTGATC"
        changed = list(query)
        changed[15] = "A" if changed[15] != "A" else "C"
        reference = "TTTT" + "".join(changed) + "AAAA"
        hits, truncated = reference_supported_alignments(query, {"ref": reference})
        self.assertFalse(truncated)
        self.assertTrue(hits)
        self.assertLessEqual(hits[0]["edit_distance"], 2)
        self.assertGreaterEqual(hits[0]["identity"], 0.90)
        reverse = query.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
        reverse_hits, _ = reference_supported_alignments(reverse, {"ref": query})
        self.assertTrue(any(hit["strand"] == "-" for hit in reverse_hits))

    def test_seeded_sample_ladder_is_nested_and_input_order_independent(self):
        read_ids = [f"read-{index}" for index in range(100)]
        subset_50 = set(sorted(read_ids, key=lambda read_id: _sample_rank(1, read_id))[:50])
        subset_100 = set(sorted(reversed(read_ids), key=lambda read_id: _sample_rank(1, read_id))[:100])
        self.assertEqual(len(subset_50), 50)
        self.assertEqual(len(subset_100), 100)
        self.assertTrue(subset_50.issubset(subset_100))

    def test_offline_replay_is_byte_deterministic(self):
        def replay_report():
            twin = AssemblyGraphTwin(k=3)
            first_batch = [Read("r1", "AAACCTG"), Read("r2", "AAACCTG")]
            second_batch = first_batch + [Read("r3", "AAATCTG")]
            twin.update(first_batch, "batch 1", source="batch-1", timestamp_utc="2000-01-01T00:00:01+00:00")
            twin.update(second_batch, "batch 2", source="batch-2", timestamp_utc="2000-01-01T00:00:02+00:00")
            return json.dumps({"snapshots": twin.snapshots}, sort_keys=True)

        first = replay_report()
        second = replay_report()
        self.assertEqual(first, second)
        report = json.loads(first)
        self.assertEqual([item["batch_read_count"] for item in report["snapshots"]], [2, 1])

    def test_optional_model_reports_missing_pytorch_clearly(self):
        if importlib.util.find_spec("torch") is not None:
            self.skipTest("PyTorch is installed; missing-dependency path does not apply")
        script = PROJECT_ROOT / "path_model.py"
        command = subprocess.run(
            [sys.executable, str(script), "train", "--train-reports", "unused.json", "--validation-reports", "unused-validation.json", "--output", "unused-model.json"],
            capture_output=True, text=True, check=False,
        )
        self.assertEqual(command.returncode, 2)
        self.assertIn("PyTorch is optional", command.stderr)


if __name__ == "__main__":
    unittest.main()

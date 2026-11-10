"""panagram.view: the app built by view() (layout, parameters, Dash endpoints) and the `panagram view` options.

The server is never started: Dash.run is replaced by a recorder, and the Flask test client
serves the layout, the callback table and the callbacks themselves.
"""

import json

import pytest

from helpers_view import components, decode_array, load_view, patch_umaps, umaps, v, view_ns  # noqa: F401

REQUIRED_IDS = {
	"genome_select_dropdown",
	"chr_select_dropdown",
	"Chrs_Info",
	"pre_regions_dropdown",
	"chr_name",
	"regional_genes",
	"step_size_out",
	"tabs",
	"all_genomes",
	"all_genomes_kmer_comp",
	"all_genome_sizes",
	"average_kmer_content",
	"pangenome_num_seqs",
	"genome_umap",
	"genome_umap_hist",
	"all_chromosomes",
	"all_chromosomes_hists",
	"Umap",
	"Umap_genome_hist",
	"chromosome",
	"primary",
	"Genes",
	"Third",
	"selected-chrom-state",
	"selected-anchor-state",
	"start-coord-state",
	"end-coord-state",
	"chr-genes-state",
}


def _by_id(app):
	return {c.id: c for c in components(app.layout) if isinstance(getattr(c, "id", None), str)}


def _texts(app):
	out = []
	for c in components(app.layout):
		child = getattr(c, "children", None)
		if isinstance(child, str):
			out.append(child)
	return out


@pytest.fixture
def open_view(anchored_index, monkeypatch):
	"""Factory: run view() with the given parameters; closes every index it opened."""
	opened = []

	def make(**kwargs):
		ns = load_view(anchored_index, monkeypatch, **kwargs)
		opened.append(ns)
		return ns

	yield make
	for ns in opened:
		ns.index.close()


# ---------------------------------------------------------------- layout


def test_layout_contains_every_component_the_callbacks_use(v):
	assert REQUIRED_IDS <= set(_by_id(v.app))


def test_header_texts_describe_the_index(v):
	texts = _texts(v.app)
	assert "Panagram of 5 genomes" in texts
	assert "K-mer length: 15" in texts
	assert "Step size: 100" in texts
	assert "Panagram" in texts


def test_genome_dropdown_lists_anchor_genomes(v):
	by_id = _by_id(v.app)
	dd = by_id["genome_select_dropdown"]
	assert [o["value"] for o in dd.options] == ["alpha", "beta", "gamma", "delta", "epsilon"]
	assert all(o["label"] == o["value"] for o in dd.options)
	assert dd.value == "alpha"
	assert dd.clearable is False


def test_chromosome_dropdown_lists_anchor_chromosomes(v):
	dd = _by_id(v.app)["chr_select_dropdown"]
	assert [o["value"] for o in dd.options] == ["chr1", "chr2"]
	assert dd.value == ""


def test_tabs_order_and_content(v):
	tabs = _by_id(v.app)["tabs"]
	assert tabs.value == "pangenome"
	assert [(t.value, t.label) for t in tabs.children] == [
		("pangenome", "Pangenome"),
		("anchor", "Anchor genome"),
		("chromosome", "Chromosome"),
	]
	ids_in = {t.value: {c.id for c in components(t.children) if isinstance(getattr(c, "id", None), str)} for t in tabs.children}
	assert ids_in["pangenome"] == {
		"all_genomes",
		"all_genomes_kmer_comp",
		"all_genome_sizes",
		"average_kmer_content",
		"pangenome_num_seqs",
	}
	assert {"genome_umap", "genome_umap_hist", "all_chromosomes", "all_chromosomes_hists"} == ids_in["anchor"]
	assert {"Umap", "Umap_genome_hist", "chromosome", "primary", "Genes", "Third"} == ids_in["chromosome"]


def test_anchor_tab_graph_height_scales_with_chromosome_count(v):
	graph = _by_id(v.app)["all_chromosomes"]
	assert graph.style["height"] == 2 * 250  # alpha has two chromosomes


def test_whole_genome_histogram_figure_is_prebuilt(v):
	"""all_chromosomes_hists is filled at start-up with the per-chromosome conservation figure."""
	graph = _by_id(v.app)["all_chromosomes_hists"]
	assert len(graph.figure.data) > 0
	assert graph.figure.layout.title.text == "K-mer conservation per chromosome"


def test_initial_state_defaults(v):
	by_id = _by_id(v.app)
	assert by_id["selected-anchor-state"].children == "alpha"
	assert by_id["selected-chrom-state"].children == "chr1"
	assert by_id["start-coord-state"].children == 0
	assert by_id["end-coord-state"].children == 2636
	assert by_id["chr_name"].children == "alpha.chr1:"


def test_graph_config_exports_svg_without_scroll_zoom(v):
	graph = _by_id(v.app)["primary"]
	assert graph.config["toImageButtonOptions"]["format"] == "svg"
	assert graph.config["scrollZoom"] is False


def test_b64_image_encodes_file_as_png_data_uri(v, tmp_path):
	path = tmp_path / "x.png"
	path.write_bytes(b"abc\x00\xff")
	assert v.b64_image(str(path)) == "data:image/png;base64,YWJjAP8="


# ---------------------------------------------------------------- parameters


def test_genome_start_end_parameters(open_view):
	ns = open_view(genome="gamma", start=100, end=900)
	by_id = _by_id(ns.app)
	assert by_id["genome_select_dropdown"].value == "gamma"
	assert by_id["selected-anchor-state"].children == "gamma"
	assert by_id["start-coord-state"].children == 100
	assert by_id["end-coord-state"].children == 900
	assert [o["value"] for o in by_id["chr_select_dropdown"].options] == ["chr1", "chr2"]


def test_default_end_is_size_of_first_chromosome(open_view):
	ns = open_view(genome="beta")
	by_id = _by_id(ns.app)
	assert by_id["start-coord-state"].children == 0
	assert by_id["end-coord-state"].children == 2886


def test_run_receives_host_port_and_debug(open_view):
	default = open_view()
	assert default.run_kwargs == {"host": "127.0.0.1", "port": "8050", "debug": False}
	custom = open_view(host="0.0.0.0", port="9123", ndebug=False)
	assert custom.run_kwargs == {"host": "0.0.0.0", "port": "9123", "debug": True}


def test_url_base_prefixes_the_endpoints(open_view):
	ns = open_view(url_base="/pan/")
	client = ns.app.server.test_client()
	assert client.get("/pan/_dash-layout").status_code == 200
	assert client.get("/pan/_dash-dependencies").status_code == 200
	assert client.get("/_dash-layout").status_code == 404


def test_max_chr_bins_parameter_reaches_the_bin_size(open_view, monkeypatch):
	ns = open_view(max_chr_bins=100)
	patch_umaps(ns, monkeypatch)
	_, main, *_ = ns.update_all_figs(1, "chr1", "alpha", 0, 0, 2636, 1)
	bar = [t for t in main.data if t.type == "bar" and t.name == "1"][0]
	assert list(bar.x)[:3] == [0, 27, 54]


def test_order_parameter_reaches_the_heatmap(open_view, monkeypatch):
	order = ["beta", "alpha", "gamma", "delta", "epsilon"]
	ns = open_view(order=order)
	patch_umaps(ns, monkeypatch)
	_, main, *_ = ns.update_all_figs(1, "chr1", "alpha", 0, 0, 2636, 1)
	heat = [t for t in main.data if t.type == "heatmap" and t.coloraxis == "coloraxis2"][0]
	assert list(heat.y) == order


def test_command_line_options_drive_view(anchored_index, monkeypatch, umaps):
	"""`panagram view <dir> gamma chr2 10 600 --ndebug ...` parsed by panagram.__main__ and run through View.run()."""
	import dash
	from simple_parsing import ArgumentParser

	from panagram.__main__ import Main

	parser = ArgumentParser(add_config_path_arg=True)
	parser.add_arguments(Main, dest="main")
	argv = [
		"view", str(anchored_index), "gamma", "chr2", "10", "600",
		"--port", "9001", "--host", "0.0.0.0", "--url_base", "/cli/", "--max_chr_bins", "50",
		"--order", "alpha", "beta", "gamma", "delta", "epsilon",
	]
	args, extra = parser.parse_known_args(argv)
	assert extra == []
	seen = {}

	def fake_run(self, **kwargs):
		seen["kwargs"] = kwargs
		seen["app"] = self

	monkeypatch.setattr(dash.Dash, "run", fake_run)
	args.main.run()
	assert seen["kwargs"] == {"host": "0.0.0.0", "port": "9001", "debug": True}
	assert seen["app"].config.url_base_pathname == "/cli/"
	by_id = _by_id(seen["app"])
	assert by_id["start-coord-state"].children == 10
	assert by_id["end-coord-state"].children == 600
	assert by_id["genome_select_dropdown"].value == "gamma"


# ---------------------------------------------------------------- Dash endpoints


def test_layout_endpoint_serves_json_with_every_id(v):
	client = v.app.server.test_client()
	resp = client.get("/_dash-layout")
	assert resp.status_code == 200
	ids = set()

	def walk(node):
		if isinstance(node, dict):
			props = node.get("props", {})
			if isinstance(props.get("id"), str):
				ids.add(props["id"])
			for value in props.values():
				walk(value)
		elif isinstance(node, list):
			for n in node:
				walk(n)

	walk(json.loads(resp.data))
	assert REQUIRED_IDS <= ids


def test_dependencies_reference_only_existing_components(v):
	"""Every input, state and output of the navigation, pangenome, anchor and chromosome callbacks is in the layout."""
	client = v.app.server.test_client()
	deps = json.loads(client.get("/_dash-dependencies").data)
	assert len(deps) == 5
	layout_ids = set(_by_id(v.app))
	checked = 0
	for dep in deps:
		if dep["output"].startswith("annotation_conservation"):
			continue  # annotation tab is disabled: neither its graph nor its radio items exist
		for item in dep["inputs"] + dep["state"]:
			assert item["id"] in layout_ids, item
		for out in dep["output"].strip(".").split("..."):
			assert out.split(".")[0] in layout_ids, out
		checked += 1
	assert checked == 4


def _post(client, dep, inputs, state=(), changed=None):
	outputs = [{"id": o.split(".")[0], "property": o.split(".")[1]} for o in dep["output"].strip(".").split("...")]
	body = {
		"output": dep["output"],
		"outputs": outputs if len(outputs) > 1 else outputs[0],
		"inputs": [{"id": i, "property": p, "value": val} for i, p, val in inputs],
		"state": [{"id": i, "property": p, "value": val} for i, p, val in state],
		"changedPropIds": [changed or f"{inputs[0][0]}.{inputs[0][1]}"],
	}
	return client.post("/_dash-update-component", json=body)


def _dep(deps, first_output):
	return [d for d in deps if d["output"].strip(".").startswith(first_output)][0]


def test_http_navigation_callback_changes_genome(v):
	client = v.app.server.test_client()
	deps = json.loads(client.get("/_dash-dependencies").data)
	dep = _dep(deps, "selected-anchor-state")
	inputs = [
		("genome_select_dropdown", "value", "gamma"),
		("chr_select_dropdown", "value", None),
		("all_chromosomes", "relayoutData", None),
		("primary", "clickData", None),
		("primary", "relayoutData", None),
		("chromosome", "selectedData", None),
		("Genes", "clickData", None),
		("Chrs_Info", "value", None),
		("pre_regions_dropdown", "value", None),
	]
	state = [
		("selected-anchor-state", "children", "alpha"),
		("selected-chrom-state", "children", "chr1"),
		("start-coord-state", "children", 0),
		("end-coord-state", "children", 2636),
		("chr-genes-state", "children", None),
	]
	resp = _post(client, dep, inputs, state)
	assert resp.status_code == 200
	body = json.loads(resp.data)["response"]
	assert body["selected-anchor-state"] == {"children": "gamma"}
	assert body["selected-chrom-state"] == {"children": "chr1"}
	assert body["start-coord-state"] == {"children": 0}
	assert body["end-coord-state"] == {"children": 2486}
	assert body["Chrs_Info"] == {"value": "0-2486"}
	assert body["chr_name"] == {"children": "gamma.chr1:"}


def test_http_pangenome_callback_returns_serialised_figures(v):
	client = v.app.server.test_client()
	deps = json.loads(client.get("/_dash-dependencies").data)
	dep = _dep(deps, "all_genomes.figure")
	resp = _post(client, dep, [("tabs", "value", "pangenome"), ("selected-anchor-state", "children", "alpha")])
	assert resp.status_code == 200
	body = json.loads(resp.data)["response"]
	assert set(body) == {"all_genomes", "all_genomes_kmer_comp", "pangenome_num_seqs", "average_kmer_content", "all_genome_sizes"}
	counts = body["pangenome_num_seqs"]["figure"]["data"][0]
	assert dict(zip(decode_array(counts["x"]), decode_array(counts["y"]))) == {"alpha": 2, "beta": 1, "gamma": 2, "delta": 1, "epsilon": 1}


def test_http_chromosome_callback_serialises_main_figure(v, umaps, synthetic_genomes):
	client = v.app.server.test_client()
	deps = json.loads(client.get("/_dash-dependencies").data)
	dep = _dep(deps, "chromosome.figure")
	inputs = [("tabs", "value", "chromosome"), ("start-coord-state", "children", 0), ("end-coord-state", "children", 2636)]
	state = [("selected-chrom-state", "children", "chr1"), ("selected-anchor-state", "children", "alpha")]
	resp = _post(client, dep, inputs, state)
	assert resp.status_code == 200
	body = json.loads(resp.data)["response"]
	assert set(body) == {"chromosome", "primary", "Third", "Genes", "Umap", "Umap_genome_hist", "regional_genes"}
	assert body["regional_genes"]["children"].startswith("Genes: alpha_chr1_private")
	bars = [t for t in body["primary"]["figure"]["data"] if t["type"] == "bar" and t.get("name") == "1"]
	assert len(bars) == 1 and (
		sum(decode_array(bars[0]["y"])) == (synthetic_genomes.oracle_bits("alpha", "chr1").sum(axis=1) == 1).sum()
	)  # k-mers private to alpha on chr1

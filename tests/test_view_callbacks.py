"""panagram.view: the four Dash callbacks (navigation, pangenome, anchor, chromosome) and the annotation helpers.

Callbacks are called as plain functions with the Dash context set to the triggering input.
"""

import dash
import numpy as np
import pandas as pd
import pytest

from helpers_view import UMAP_ROWS, call, umaps, v, view_ns  # noqa: F401

NO = dash.no_update


def nav(v, trigger, anchor="alpha", chrom="chr1", start=0, end=2636, **inputs):
	"""Call nav_callback; `inputs` are the callback arguments by keyword (default None)."""
	names = [
		"anchor_dropdown",
		"chr_dropdown",
		"anctab_chrs_relayout",
		"chrtab_primary_click",
		"chrtab_primary_relayout",
		"chrtab_chr_select",
		"chrtab_gene_click",
		"user_chr_coords",
		"pre_selected_region",
	]
	args = [inputs.pop(n, None) for n in names]
	assert not inputs, inputs
	return call(v.nav_callback, trigger, *args, anchor, chrom, start, end, None)


# outputs: anchor, chrom, start, end, chromosome.relayoutData, Chrs_Info.value, chr_name.children


def test_nav_genome_dropdown_resets_to_first_chromosome(v):
	out = nav(v, "genome_select_dropdown.value", anchor_dropdown="gamma")
	assert out == ("gamma", "chr1", 0, 2486, None, "0-2486", "gamma.chr1:")


def test_nav_genome_dropdown_to_genome_with_other_size(v):
	out = nav(v, "genome_select_dropdown.value", anchor_dropdown="beta")
	assert out == ("beta", "chr1", 0, 2886, None, "0-2886", "beta.chr1:")


def test_nav_chromosome_dropdown(v):
	out = nav(v, "chr_select_dropdown.value", chr_dropdown="chr2")
	assert out == ("alpha", "chr2", 0, 1686, None, "0-1686", "alpha.chr2:")


def test_nav_chromosome_dropdown_unknown_name_gives_first(v):
	out = nav(v, "chr_select_dropdown.value", chr_dropdown="chrZ")
	assert out == ("alpha", "chr1", 0, 2636, None, "0-2636", "alpha.chr1:")


def test_nav_pre_selected_region(v):
	out = nav(v, "pre_regions_dropdown.value", pre_selected_region="chr2:100-900")
	assert out == ("alpha", "chr2", 100, 900, None, "100-900", "alpha.chr2:")


def test_nav_pre_selected_region_out_of_range_is_ignored(v):
	"""End beyond the chromosome: coordinates are not updated and the previous region is shown."""
	out = nav(v, "pre_regions_dropdown.value", pre_selected_region="chr1:100-99999")
	assert out[0] == "alpha" and out[1] == "chr1"
	assert out[2] is NO and out[3] is NO
	assert out[5:] == ("0-2636", "alpha.chr1:")


def test_nav_pre_selected_region_none_changes_nothing(v):
	out = nav(v, "pre_regions_dropdown.value", pre_selected_region=None)
	assert out[0] == "alpha" and out[1] == "chr1"
	assert out[2] is NO and out[3] is NO
	assert out[5] == "0-2636"


def test_nav_typed_coordinates(v):
	out = nav(v, "Chrs_Info.value", user_chr_coords="100-500")
	assert out == ("alpha", "chr1", 100, 500, None, "100-500", "alpha.chr1:")
	padded = nav(v, "Chrs_Info.value", user_chr_coords=" 200-300 ")
	assert padded[2:4] == (200, 300)


def test_nav_typed_coordinates_end_at_chromosome_size_is_accepted(v):
	"""chr1 of alpha is 2636 long: an end of 2636 is valid, 2637 is not."""
	ok = nav(v, "Chrs_Info.value", user_chr_coords="0-2636")
	assert ok[2:4] == (0, 2636)
	too_long = nav(v, "Chrs_Info.value", user_chr_coords="0-2637")
	assert too_long[2] is NO and too_long[3] is NO


@pytest.mark.parametrize("typed", ["500-100", "100-99999"])
def test_nav_typed_coordinates_invalid_range_keeps_previous(v, typed):
	out = nav(v, "Chrs_Info.value", user_chr_coords=typed)
	assert out[2] is NO and out[3] is NO
	assert out[5] == "0-2636"
	assert out[1] == "chr1"


def test_nav_gene_click_sets_gene_bounds(v):
	click = {"points": [{"customdata": [301, 900], "text": "alpha_chr1_CORE_1"}]}
	out = nav(v, "Genes.clickData", chrtab_gene_click=click)
	assert out[2:4] == (301, 900)
	assert out[5] == "301-900"


def test_nav_gene_click_none_changes_nothing(v):
	out = nav(v, "Genes.clickData", chrtab_gene_click=None)
	assert out[2] is NO and out[3] is NO and out[5] == "0-2636"


@pytest.mark.parametrize("axis", ["x", "x2", "x3"])
def test_nav_chromosome_selection_ranges(v, axis):
	selected = {"range": {axis: [120.7, 480.2]}}
	out = nav(v, "chromosome.selectedData", chrtab_chr_select=selected)
	assert out[2:4] == (120, 480)
	assert out[4] == selected
	assert out[5] == "120-480"


def test_nav_primary_relayout_zoom(v):
	zoom = {"xaxis2.range[0]": 300.4, "xaxis2.range[1]": 700.9}
	out = nav(v, "primary.relayoutData", chrtab_primary_relayout=zoom)
	assert out[2:4] == (300, 700)
	assert out[5] == "300-700"


def test_nav_primary_relayout_other_events_change_nothing(v):
	out = nav(v, "primary.relayoutData", chrtab_primary_relayout={"autosize": True})
	assert out[2] is NO and out[3] is NO and out[5] == "0-2636"
	none = nav(v, "primary.relayoutData", chrtab_primary_relayout=None)
	assert none[2] is NO and none[3] is NO


def test_nav_anchor_tab_event_without_axis_keys_changes_nothing(v):
	out = nav(v, "all_chromosomes.relayoutData", anctab_chrs_relayout={"autosize": True})
	assert out[0] is NO
	assert out[1] == "chr1" and out[2] is NO and out[3] is NO
	assert out[5] == "0-2636"


def test_nav_without_trigger_keeps_state(v):
	"""The initial call (no triggering input) leaves coordinates and genome as they are."""
	out = call(v.nav_callback, None, *([None] * 9), "gamma", "chr2", 0, 686, None)
	assert out[0] == "gamma" and out[1] == "chr2"
	assert out[2] is NO and out[3] is NO
	assert out[5:] == ("0-686", "gamma.chr2:")


def test_nav_output_count_matches_declared_outputs(v):
	keys = [k for k in v.app.callback_map if "selected-anchor-state.children" in k.split("...")[0]]
	assert len(keys) == 1
	n_outputs = len(keys[0].strip(".").split("..."))
	out = nav(v, "genome_select_dropdown.value", anchor_dropdown="alpha")
	assert len(out) == n_outputs == 7


# ---------------------------------------------------------------- tab callbacks


def _trace_xy(fig):
	return list(fig.data[0].x), list(fig.data[0].y)


def test_pangenome_callback_other_tab_returns_empty_figures(v):
	out = call(v.pangenome_callback, "tabs.value", "chromosome", "alpha")
	assert out == ({},) * 5


@pytest.mark.parametrize("trigger", ["tabs.value", "selected-anchor-state.children"])
def test_pangenome_callback_builds_five_figures(v, trigger):
	dend, comp, count, avgs, sizes = call(v.pangenome_callback, trigger, "pangenome", "delta")
	assert dend is v.all_genomes_dend_fig
	assert comp is v.pangenome_comp_fig
	assert _trace_xy(count) == _trace_xy(v.make_genome_count_plot("delta"))
	assert dict(zip(*_trace_xy(count))) == {"alpha": 2, "beta": 1, "gamma": 2, "delta": 1, "epsilon": 1}
	assert _trace_xy(avgs) == _trace_xy(v.plot_avgs("delta"))
	assert _trace_xy(sizes) == _trace_xy(v.make_genome_size_plot("delta"))
	assert [s.x0 for s in count.layout.shapes] == ["delta"]


def test_pangenome_callback_unrelated_trigger_is_no_update(v):
	out = call(v.pangenome_callback, "start-coord-state.children", "pangenome", "alpha")
	assert all(o is NO for o in out) and len(out) == 5


def test_anchor_callback_other_tab_gives_empty_figures_and_options(v):
	*figs, options = call(v.anchor_callback, "tabs.value", "pangenome", "gamma")
	assert figs == [{}, {}, {}]
	assert options == [{"label": "chr1", "value": "chr1"}, {"label": "chr2", "value": "chr2"}]


def test_anchor_callback_anchor_tab(v):
	whole, umap, hist, options = call(v.anchor_callback, "tabs.value", "anchor", "beta")
	heat = [t for t in whole.data if t.type == "heatmap"]
	assert len(heat) == 3  # one chromosome: cluster, universal, private rows
	assert whole.layout.height == 250
	assert len(umap.data[0].x) == len(v.index.genomes["beta"].genome_umap)
	assert list(hist.data[0].y) == [0]
	assert options == [{"label": "chr1", "value": "chr1"}]


def test_anchor_callback_anchor_change_updates_chromosome_options(v):
	for anchor, chroms in (("alpha", ["chr1", "chr2"]), ("epsilon", ["chr1"])):
		options = call(v.anchor_callback, "selected-anchor-state.children", "anchor", anchor)[3]
		assert [o["value"] for o in options] == chroms


def test_chromosome_callback_other_tab_is_no_update(v):
	out = call(v.chromosome_callback, "tabs.value", "anchor", 0, 2636, "chr1", "alpha")
	assert len(out) == 7 and all(o is NO for o in out)


def test_chromosome_callback_matches_update_all_figs(v, umaps):
	out = call(v.chromosome_callback, "tabs.value", "chromosome", 100, 900, "chr1", "alpha")
	ref = v.update_all_figs(1, "chr1", "alpha", 0, 100, 900, v.index.lowres_step)
	assert len(out) == 7
	assert out[6] == ref[6]
	assert out[6].rstrip(", ") == "Genes: alpha_chr1_private, alpha_chr1_CORE_1"
	for got, want in zip(out[:6], ref[:6]):
		assert [t.type for t in got.data] == [t.type for t in want.data]
		assert [getattr(t, "name", None) for t in got.data] == [getattr(t, "name", None) for t in want.data]
	assert out[0].layout.xaxis.range == (0, 2636)  # the whole-chromosome overview
	assert out[1].layout.xaxis2.range == (100, 900)  # the zoomed main figure


def test_chromosome_callback_region_drives_gene_text(v, umaps):
	one = call(v.chromosome_callback, "tabs.value", "chromosome", 1200, 1500, "chr1", "alpha")
	assert one[6] == "Genes: alpha_chr1_PAIR_3"
	none = call(v.chromosome_callback, "tabs.value", "chromosome", 2300, 2600, "chr1", "alpha")
	assert none[6] == ""


# ---------------------------------------------------------------- annotation helpers


def _annotation_file(tmp_path):
	rows = [
		("geneA", "chr1", 100, 400, 2.5, 0.5, "ID=geneA"),
		("geneB", "chr2", 1000, 1100, 4.0, 1.5, "ID=geneB;Name=B"),
		("geneC", "chr1", 50, 450, 1.0, 0.0, "ID=geneC"),
	]
	path = tmp_path / "annotation.tsv"
	path.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))
	return path


def test_annotation_tab_info_reads_columns(v, tmp_path):
	df = v.annotation_tab_info(str(_annotation_file(tmp_path)))
	assert list(df["Name"]) == ["geneA", "geneB", "geneC"]
	assert list(df["chr"]) == ["chr1", "chr2", "chr1"]
	assert list(df["start"]) == [100, 1000, 50]
	assert list(df["end"]) == [400, 1100, 450]
	assert list(df["Length"]) == [300, 100, 400]
	assert list(df["Average"]) == [2.5, 4.0, 1.0]
	assert list(df["Standard_dev"]) == [0.5, 1.5, 0.0]
	assert list(df["attr"]) == ["ID=geneA", "ID=geneB;Name=B", "ID=geneC"]


def test_annotation_tab_plot_log2_length_colour(v, tmp_path):
	df = v.annotation_tab_info(str(_annotation_file(tmp_path)))
	fig = v.annotation_tab_plot(df, "Length", 1)
	main = fig.data[0]
	assert list(main.x) == [2.5, 4.0, 1.0]
	assert list(main.y) == [0.5, 1.5, 0.0]
	assert list(main.marker.color) == pytest.approx(list(np.log2([300, 100, 400])))
	assert list(main.hovertext) == ["geneA", "geneB", "geneC"]
	assert fig.layout.clickmode == "event+select"


def test_annotation_tab_plot_plain_length_colour(v, tmp_path):
	df = v.annotation_tab_info(str(_annotation_file(tmp_path)))
	fig = v.annotation_tab_plot(df, "Length", 0)
	assert list(fig.data[0].marker.color) == [300, 100, 400]


def test_annotation_tab_plot_by_chromosome_groups_points(v, tmp_path):
	df = v.annotation_tab_info(str(_annotation_file(tmp_path)))
	fig = v.annotation_tab_plot(df, "chr", 0)
	groups = {t.name: list(t.hovertext) for t in fig.data if t.type == "scatter" and t.hovertext is not None}
	assert groups == {"chr1": ["geneA", "geneC"], "chr2": ["geneB"]}


@pytest.mark.parametrize("radio", [" Log2 gene length", " Gene length", " Chromosome", "Bed file", "other"])
def test_annotation_tab_callback_returns_blank_figure(v, radio):
	"""pins the current output: the callback computes the colour mode but returns an empty figure (the tab is disabled)."""
	fig = call(v.annotation_tab_callback, "radios.value", radio)
	assert len(fig.data) == 0

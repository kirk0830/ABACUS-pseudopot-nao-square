"""
Driver Name: EcutwfcConv_20240319

Driver description:
This is driver for APNS pages pseudopotential convergence test
Two kinds of data are post-processed:
1. ecutwfc - energy per atom
2. ecutwfc - pressure of whole cell
"""

import argparse
def entry():
    parser = argparse.ArgumentParser(description="APNS pseudopotential convergence test")
    # add -i
    parser.add_argument("-i", "--input", type=str, help="input json file")
    args = parser.parse_args()

    return args.input

import json
def testname_pspotid(testname: str):
    """return pspotid from testname"""
    def translate_kind(testname):
        """translate pspotid to real pseudopotential name"""
        dictionary = {"dojo": "PseudoDojo", "pslnc": "PSlibrary (Norm-Conserving)", 
                      "pslrrkjus": "PSlibrary (RRKJUS)", "pslncpaw": "PSlibrary (PAW)",
                      "GTH": "Goedecker-Teter-Hutter", "HGH": "Hartwigsen-Goedecker-Hutter"}
        return dictionary[testname] if testname in dictionary.keys() else testname.upper()
    def translate_version(kind, version):
        if kind.upper() == "GTH":
            return "-"+version
        elif kind.upper() == "PD":
            return version
        else:
            version = ".".join([str(version[i]) for i in range(len(version))])
            return "-v" + version if version != "" else ""
    # TEMPORARY FIX
    # this is because presently only GTH LnPP1 are collected, while there are indeed other versions
    # of GTH pseudopotentials, so we need to fix the name temporarily. Other versions like UZH, LnPP2,
    # ... will be added in the future.
    testname = "gthLnPP1" if testname == "gth" else testname
    with open("download/pseudopotentials/description.json", "r") as f:
        description = json.load(f)
    for key in description.keys():
        if key.replace("_", "") == testname:
            val = description[key]
            with open(val + "/description.json", "r") as f:
                data = json.load(f)
            data = {key: data[key] for key in data.keys() if key != "files"}
            kind = translate_kind(data["kind"])
            version = translate_version(data["kind"], data["version"])
            label = kind + version
            label += " (" + data["appendix"] + ")" if data["appendix"] != "" else ""
            label = label.strip()
            return label, key
    print("WARNING: testname", testname, "not found in description.json")
    return None, None

import apns.module_pseudo.parse as ampp
def z_valence(element: str, pspotid: str):
    """pspotid here should be the "key" returned by function testname_pspotid"""
    with open("download/pseudopotentials/description.json", "r") as f:
        description = json.load(f)
    path = description[pspotid]
    path += "/" if path[-1] != "/" else ""
    fdes = path + "description.json"
    with open(fdes, "r") as f:
        description = json.load(f)
    fpspot = path + description["files"][element]
    return ampp.z_valence(fpspot)

def categorize_byelement(labels, data):
    """the data will arrange like:
    ```python
    (
        [('Ag', '8566', 'dojo05', 90.0), ('Ag', '8566', 'pd03', 90.0), 
        ('Ag', '8566', 'sg1512', 60.0), ('Ag', '8566', 'sg1510', 60.0), 
        ...
        ],
        [
            [
                [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 150.0, 200.0, 300.0], 
                [111.00118765492425, 33.39993995924942, 9.262190618749628, 2.1097572079243037, 0.3720193069993911, 0.041972779749812617, 
                0.002425373099868011, 0.0005534015244847978, 0.00044674614946416114, 0.00019317457463330356, 4.175727463007206e-05, 0.0]
            ], 
            [
                [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 150.0, 200.0, 300.0], 
                [111.14227980914984, 33.32067861827545, 9.212816671300061, 2.113927723249617, 0.3799892410997927, 0.04445111024961079, 
                0.0027109032753287465, 0.000579334774556628, 0.00048263027474604314, 0.00022845432522444753, 4.262575021130033e-05, 0.0]
            ], 
            ...
        ],
        [('Ag', '8566', 'dojo05', 90.0), ('Ag', '8566', 'pd03', 90.0), 
        ('Ag', '8566', 'sg1512', 60.0), ('Ag', '8566', 'sg1510', 60.0), 
        ...
        ],
        [
            [
                [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 150.0, 200.0, 300.0], 
                [111.00118765492425, 33.39993995924942, 9.262190618749628, 2.1097572079243037, 0.3720193069993911, 0.041972779749812617, 
                0.002425373099868011, 0.0005534015244847978, 0.00044674614946416114, 0.00019317457463330356, 4.175727463007206e-05, 0.0]
            ], 
            [
                [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 150.0, 200.0, 300.0], 
                [111.14227980914984, 33.32067861827545, 9.212816671300061, 2.113927723249617, 0.3799892410997927, 0.04445111024961079, 
                0.0027109032753287465, 0.000579334774556628, 0.00048263027474604314, 0.00022845432522444753, 4.262575021130033e-05, 0.0]
            ], 
            ...
        ]
    )
    ```
    each time import two elements, the first is labels, the second is data.

    This function will categorize the data by element, 
    pspotids first indiced by element, second by pspotid,
    the same for decomposed_x and decomposed_y

    Return:
    elements: list of elements
    pspotids: list of list of pseudopotential ids
    decomposed_x: list of list of x data
    decomposed_y: list of list of y data
    """
    assert len(labels) == len(data)
    for label in labels:
        assert len(label) == 4 # element, mpid, pspotid, ecutwfc
        element, mpid, pspotid, ecutwfc = label
        assert isinstance(element, str)
        assert isinstance(mpid, str)
        assert isinstance(pspotid, str)
        assert isinstance(ecutwfc, float)
    for piece in data:
        assert len(piece) == 2 # ecutwfc, data
        assert isinstance(piece[0], list)
        assert isinstance(piece[1], list)
        assert len(piece[0]) == len(piece[1])

    elements = [] # is designed to only collect unique elements
    pspotids = [] # for each element, it would be a list of pspotids
    decomposed_x = [] # corresponds to pspotids, for each element, there would be several tests, and each test has a series of ecutwfc
    decomposed_y = [] # corresponds to pspotids, for each element, there would be several tests, and each test has a series of data
    for itst in range(len(labels)): # for each test, loop. one test is defined as fixed pseudopotential, series of SCF at different ecutwfc
        # there always are duplicates for element, because four keys define the uniqueness
        element, mpid, pspotid, ecutwfc = labels[itst] # get element, mpid, pspotid, and converged ecutwfc information
        if element not in elements: # if element is not in elements, append it
            elements.append(element)
            pspotids.append([]) # store a new element, thus append a new list for it
            decomposed_x.append([]) # store a new element, thus append a new list for it
            decomposed_y.append([]) # store a new element, thus append a new list for it
        idx = elements.index(element)
        pspotids[idx].append(pspotid)
        decomposed_x[idx].append(data[itst][0])
        decomposed_y[idx].append(data[itst][1])

    return elements, pspotids, decomposed_x, decomposed_y

import apns.module_analysis.postprocess.conv.ecutwfc_eks as amapeks
import apns.module_analysis.postprocess.conv.ecutwfc_press as amapprs
def search(path: str, ethr: float = 1e-3, pthr: float = 0.1):

    return amapeks.run(path, ethr), amapprs.run(path, pthr)

def run():
    # always make independent on path for run() function, instead, use entry()
    path = entry()
    # postprocess
    (conv_eks, data_eks), (conv_prs, data_prs) = search(path)
    # categorize by element
    elements, pspotids, eks_x, eks_y = categorize_byelement(conv_eks, data_eks)
    elements, pspotids, prs_x, prs_y = categorize_byelement(conv_prs, data_prs)

    conv_results = {}
    for iconv in range(len(conv_eks)):
        if not conv_eks[iconv][0] in conv_results.keys():
            conv_results[conv_eks[iconv][0]] = []
        conv_results[conv_eks[iconv][0]].append((conv_eks[iconv][2], max(conv_eks[iconv][3], conv_prs[iconv][3])))

    # remember: eks_x/prs_x/eks_y/prs_y is indiced by [element][pspotid][ecutwfc] to get a float value
    for i in range(len(elements)):
        element = elements[i]
        _pspotnames, _pspotids = zip(*[testname_pspotid(id) for id in pspotids[i]])
        # DESIGN NOTE OF LOG_PLOTS AND STACK_LINEPLOTS
        # -------------------------------------------
        # LOG_PLOTS:
        # in some case, Bohrium crashes when calculating the pressure, this will lead to absent of a pair of (ecutwfc, pressure).
        # therefore, the length of eks_x and prs_x may not be the same, (the same for eks_y and prs_y).
        # In the following discrete_logplots() function, xs and ys support the input arranged in a way such that xs and ys are indiced
        # by [pspotid][ecutwfc], thus the eks_x[i][j] would be of meaning of ecutwfc list for j-th pseudopotential for element i.

        # what want to draw is, two subplots, one for energy convergence, one for pressure convergence.
        # For each subplot, there are several lines, each line corresponds to one pseudopotential. All x-y pairs are not needed
        # to be the same length, but for each pair, x and y should have identical length.
        # Therefore the most natural way to organize the data would be, first indiced by suplot, second by line and last by data.
        # the xs and ys should be instead indiced by [subplots][lines][data]
        
        # eks_x[i]/prs_x[i]/eks_y[i]/prs_y[i], indiced by [pspotid][ecutwfc]
        # or say shape = (npspotid, necutwfc)
        # let them be (nproperties, npspotid, necutwfc)
        xs = [eks_x[i], prs_x[i]]
        ys = [eks_y[i], prs_y[i]]

        logplot_style = {"highlight_ys": [1e-3, 0.1], "nrows": 1, 
                         "xtitle": "Kinetic energy cutoff for planewave basis (ecutwfc, in Ry)", 
                         "ytitle": ["Absolute Kohn-Sham energy difference per atom (eV)", 
                                    "Absolute pressure difference (kbar)"], 
                         "suptitle": element, 
                         "supcomment": "NOTE: Semitransparent lines correspond to red circles in the lineplot, \
with precision threshold 1e-3 eV/atom and 0.1 kbar respectively. Absence of data points result from SCF\
 convergence failure or walltime limit.",
                         "labels": _pspotnames, "fontsize": 12}
        fig, ax = discrete_logplots(xs=xs, ys=ys, **logplot_style)
        plt.savefig(f"{element}_logplot.svg")
        plt.close()
        # STACK_LINEPLOTS
        # stack lineplots are organized in a way that, for each element, first indiced by pseudopotential,
        # then by property, and last by data. Therefore, the xs and ys should be indiced by [pspotid][property][data]
        # the shape of xs and ys should be (npspotid, nproperties, ndata)

        # eks_x[i]/prs_x[i]/eks_y[i]/prs_y[i], indiced by [pspotid][ecutwfc]
        # or say shape = (npspotid, necutwfc)
        # let them be (npspotid, nproperties, necutwfc)
        xs = [[eks_x[i][j], prs_x[i][j]] for j in range(len(_pspotnames))]
        ys = [[eks_y[i][j], prs_y[i][j]] for j in range(len(_pspotnames))]

        lineplot_style = {"highlight_xs": conv_results[element], "ncols": 1, 
                          "subtitles": _pspotnames, 
                          "z_vals": [z_valence(element, _pspotid) for _pspotid in _pspotids], 
                          "grid": True,
                          "xtitle": "Kinetic energy cutoff for planewave basis (ecutwfc, in Ry)", 
                          "ytitle": ["Kohn-Sham energy difference per atom (eV)", 
                                     "Pressure difference (kbar)"],
                          "suptitle": element, 
                          "supcomment": "NOTE: The red circle indicates the converged ecutwfc wrt. ecutwfc$_{max}$\
 with precision threshold 1e-3 eV/atom and 0.1 kbar respectively. Absence of data points result from SCF convergence\
 failure or walltime limit.",
                          "fontsize": 12}
        fig, ax = stack_lineplots(xs=xs, ys=ys, **lineplot_style)
        plt.savefig(f"{element}.svg")
        plt.close()
        import apns.module_analysis.external_frender.htmls as amaeh
        html = amaeh.pseudopotentials(element=element, 
                                      xc_functional="PBE", 
                                      software="ABACUS",
                                      fconv=f"{element}.svg",
                                      fconvlog=f"{element}_logplot.svg")
        with open(f"{element}.md", "w") as f:
            f.write(html)

def styles_factory(property: str = "color", val: float = -1, ndim: int = None) -> list:
    if property == "color":
        colorpool = ["#2A306A", "#24B5A5", "#1DB8A8", "#015BBC", "#EACE4F"]
        return [colorpool[i % len(colorpool)] for i in range(ndim)]
    elif property == "marker":
        markerpool = ["o", "s", "D", "v", "^", "<", ">", "p", "P", "*", "h", "H", "+", "x", "X", "|", "_"]
        return [markerpool[i % len(markerpool)] for i in range(ndim)]
    elif property == "markersize":
        return [5 if val < 0 else val] * ndim
    elif property == "linestyle":
        linestylepool = ["-", "--", "-.", ":"]
        return [linestylepool[i % len(linestylepool)] for i in range(ndim)]
    elif property == "linewidth":
        return [1 if val < 0 else val] * ndim
    elif property == "alpha":
        return [1.0 if val < 0 else val] * ndim
    else:
        raise ValueError("Unknown property")

import matplotlib.pyplot as plt
import numpy as np
def stack_lineplots(xs: list, ys: list, **kwargs):
    nsubplts = len(xs)
    nlines = len(ys[0])
    for i in range(1, nsubplts):
        assert len(xs[i]) == nlines
        assert len(ys[i]) == nlines

    highlight_xs = kwargs.get("highlight_xs", None)
    ncols = kwargs.get("ncols", 1)
    subtitles = kwargs.get("subtitles", None)
    labels = kwargs.get("labels", None)
    colors = kwargs.get("colors", None)
    markers = kwargs.get("markers", None)
    markersizes = kwargs.get("markersizes", None)
    linestyles = kwargs.get("linestyles", None)
    linewidths = kwargs.get("linewidths", None)
    alpha = kwargs.get("alpha", 1.0)
    grid = kwargs.get("grid", True)
    fontsize = kwargs.get("fontsize", 12)
    z_vals = kwargs.get("z_vals", None)
    subplotsize = kwargs.get("subplotsize", (20, 1))

    npspots = nsubplts
    subtitles = ["subplot " + str(i) for i in range(npspots)] if subtitles is None else subtitles
    assert len(subtitles) == npspots
    nprptys = nlines
    labels = ["data " + str(i) for i in range(nprptys)] if labels is None else labels
    assert len(labels) == nprptys

    colors = styles_factory(property="color", ndim=nprptys) if colors is None else colors
    assert len(colors) == nprptys
    markers = styles_factory(property="marker", ndim=nprptys) if markers is None else markers
    assert len(markers) == nprptys
    markersizes = styles_factory(property="markersize", ndim=nprptys) if markersizes is None else markersizes
    assert len(markersizes) == nprptys
    linestyles = styles_factory(property="linestyle", ndim=nprptys) if linestyles is None else linestyles
    assert len(linestyles) == nprptys
    linewidths = styles_factory(property="linewidth", ndim=nprptys) if linewidths is None else linewidths
    assert len(linewidths) == nprptys
    alpha = styles_factory(property="alpha", val=alpha, ndim=nprptys)
    assert len(alpha) == nprptys

    assert z_vals is None or len(z_vals) == npspots
    # create figure and axes
    nrows = nsubplts // ncols + (nsubplts % ncols > 0)
    fig, ax = plt.subplots(nrows, ncols, figsize=(subplotsize[0] * ncols, subplotsize[1] * nrows), squeeze=False)
    
    xtitle_styles = {"ha": "center", "va": "center", "transform": fig.transFigure, "fontsize": fontsize}
    ytitle_styles = {"ha": "center", "va": "center", "rotation": "vertical", "transform": fig.transFigure, "fontsize": fontsize}

    # for there may be data failed to converge or due to failure of Bohrium platform, the length of
    # xs between different pseudopotentials may not be the same, so we need to find the minimum and maximum
    # of ecutwfc for each pseudopotential

    ecutwfc_min = min([min([min(xs[i][j]) for j in range(nprptys)]) for i in range(npspots)])
    ecutwfc_max = max([max([max(xs[i][j]) for j in range(nprptys)]) for i in range(npspots)])

    twinxs = []
    # plot
    for i in range(npspots): # loop over all pseudopotentials (although I am not willing to make it too specific)
                              # , but it is a rather easy way to explain the code
        row, col = i // ncols, i % ncols
        twinxs.append([])
        for j in range(nprptys):
            style = {"color": colors[j], "marker": markers[j], "markersize": markersizes[j],
                     "linestyle": linestyles[j], "linewidth": linewidths[j], "alpha": alpha[j]}
            twinxs[i].append(ax[row, col] if j == 0 else ax[row, col].twinx())
            twinxs[i][j].plot(xs[i][j], ys[i][j], label=labels[j], **style)
            # add yticks for each line
            ylim = np.max(np.abs(ys[i][j]))
            yticks = [-ylim*0.75, 0, ylim*0.75]
            ylims = [-ylim*1.25, ylim*1.25]
            # set yticks with color in the same color as the line
            twinxs[i][j].set_yticks(yticks)
            twinxs[i][j].set_yticklabels([f"{yticks[i]:.2f}" for i in range(len(yticks))], color=colors[j])
            # change color of y axis
            twinxs[i][j].spines["right" if j == 1 else "left"].set_color(colors[j])
            twinxs[i][j].set_ylim(ylims[0], ylims[1])
            twinxs[i][j].set_xlim(ecutwfc_min - 10, ecutwfc_max + 10)

        # add grid
        twinxs[i][0].grid(grid)
        # add subtitle at right, text left aligned
        pspotid_style = {"horizontalalignment": "right", "verticalalignment": "top", 
                         "transform": twinxs[i][0].transAxes, "fontsize": fontsize,
                         "backgroundcolor": "white"}
        twinxs[i][0].text(0.995, 0.9, subtitles[i], **pspotid_style)
        twinxs[i][0].text(0.995, 0.3, "$Z$ = %d"%z_vals[i], **pspotid_style) if z_vals is not None else None
        twinxs[i][0].axhline(0, color="black", linewidth=2, alpha=0.1)

        conv_marker_style = {"markersize": 15, "markerfacecolor": "none", "markeredgecolor": "red", "markeredgewidth": 2,
                             "zorder": 10, "alpha": 0.5}
        if highlight_xs is not None:
            # add a circle at the x position, on the topmost layer
            twinxs[i][0].plot(highlight_xs[i][1], 0, "o", **conv_marker_style)

    # subplot size adjustment, leave no space between subplots vertically
    plt.subplots_adjust(hspace=0.0, wspace=0.2)

    # xtitle, only add xtitle to the last row
    xtitle = kwargs.get("xtitle", None)
    if xtitle is not None:
        plt.text(0.5, 0.05, xtitle, **xtitle_styles)
    # ytitle, only add ytitle to the left-middle
    ytitle = kwargs.get("ytitle", None)
    if ytitle is not None:
        ytitle = [ytitle] if isinstance(ytitle, str) else ytitle
        for i in range(len(ytitle)):
            x = 0.075 if i == 0 else 0.95 + (i - 1)*0.05
            plt.text(x, 0.5, ytitle[i], **ytitle_styles)
    # suptitle
    suptitle = kwargs.get("suptitle", None)
    if suptitle is not None:
        plt.suptitle(suptitle, fontsize=fontsize*1.5)
    # supcomment, below the suptitle
    supcomment = kwargs.get("supcomment", None)
    if supcomment is not None:
        supcomment_style = {"ha": "center", "va": "center", "transform": fig.transFigure, 
                            "fontsize": fontsize, "style": "italic"}
        plt.text(0.5, 0.91, supcomment, **supcomment_style)

    # set overall fontstyle
    plt.rcParams["font.family"] = "Arial"
    return fig, ax

def discrete_logplots(xs: list, ys: list, **kwargs):

    nsubplts = len(xs)
    nlines = len(ys[0])
    for i in range(1, nsubplts):
        assert len(xs[i]) == nlines
        assert len(ys[i]) == nlines
    
    xtitle = kwargs.get("xtitle", None)
    ytitle = kwargs.get("ytitle", None)
    labels = kwargs.get("labels", None)
    colors = kwargs.get("colors", None)
    markers = kwargs.get("markers", None)
    markersizes = kwargs.get("markersizes", None)
    linestyles = kwargs.get("linestyles", None)
    linewidths = kwargs.get("linewidths", None)
    alpha = kwargs.get("alpha", 1.0)
    grid = kwargs.get("grid", True)
    fontsize = kwargs.get("fontsize", 12)
    subplotsize = kwargs.get("subplotsize", (10, 10))
    nrows = kwargs.get("nrows", 1)

    highlight_ys = kwargs.get("highlight_ys", None)

    npspots = nlines
    labels = ["data " + str(i) for i in range(npspots)] if labels is None else labels
    assert len(labels) == npspots
    nprptys = nsubplts
    ytitle = ["subplot " + str(i) for i in range(nprptys)] if ytitle is None else ytitle
    if isinstance(ytitle, str):
        ytitle = [ytitle] * nprptys
    assert len(ytitle) == nprptys

    colors = styles_factory(property="color", ndim=npspots) if colors is None else colors
    assert len(colors) == npspots
    markers = styles_factory(property="marker", ndim=npspots) if markers is None else markers
    assert len(markers) == npspots
    markersizes = styles_factory(property="markersize", ndim=npspots) if markersizes is None else markersizes
    assert len(markersizes) == npspots
    linestyles = styles_factory(property="linestyle", ndim=npspots) if linestyles is None else linestyles
    assert len(linestyles) == npspots
    linewidths = styles_factory(property="linewidth", ndim=npspots) if linewidths is None else linewidths
    assert len(linewidths) == npspots
    alpha = styles_factory(property="alpha", val=alpha, ndim=npspots)
    assert len(alpha) == npspots
    
    # create figure and axes
    ncols = nprptys // nrows + (nprptys % nrows > 0)
    fig, ax = plt.subplots(nrows, ncols, figsize=(subplotsize[0] * ncols, subplotsize[1] * nrows), squeeze=False)

    # plot
    for i in range(nprptys): # for each property
        for j in range(npspots): # for each pseudopotential
            # because the last value of y is always to be 0, so it is not plotted
            styles = {"color": colors[j], "marker": markers[j], "markersize": markersizes[j], 
                      "linestyle": linestyles[j], "linewidth": linewidths[j], "alpha": alpha[j]}
            ax[0, i].plot(xs[i][j][:-1], np.abs(np.array(ys[i][j]))[:-1], label=labels[j], **styles)
            ax[0, i].set_yscale("log")
            ax[0, i].grid(grid)
            ax[0, i].legend(fontsize=fontsize, shadow=True)
            # set xtitle
            ax[0, i].set_xlabel(xtitle, fontsize=fontsize)
            ax[0, i].set_ylabel(ytitle[i], fontsize=fontsize)
        if highlight_ys is not None:
            ax[0, i].axhline(highlight_ys[i], color="red", linewidth=5, alpha=0.1)

    # suptitle
    suptitle = kwargs.get("suptitle", None)
    if suptitle is not None:
        plt.suptitle(suptitle, fontsize=fontsize*1.5)
    # supcomment, below the suptitle
    supcomment = kwargs.get("supcomment", None)
    if supcomment is not None:
        # italicize the supcomment
        supcomment_style = {"ha": "center", "va": "center", "transform": fig.transFigure, 
                            "fontsize": fontsize, "style": "italic"}
        plt.text(0.5, 0.925, supcomment, **supcomment_style)

    # set overall fontstyle
    plt.rcParams["font.family"] = "Arial"
    return fig, ax

if __name__ == "__main__":
    # this should not be changed no matter what kind of postprocess is!
    run()
"""this file is for greping and rendering svg plots of ecutwfc convergence test for new version of APNS"""

def collect(folder: str):
    print("* * * Collect ABACUS result * * *".center(100))
    import apns.analysis.postprocess.read_abacus_out as read_abacus_out
    import apns.pspot.parse as ppparse
    import os, re
    from apns.analysis.apns2_utils import read_apnsjob_desc, convert_fpp_to_ppid
    result = {}
    for root, _, files in os.walk(folder):
        for file in files:
            if re.match(r"(running_)(\w+)(\.log)", file): # reach the inner most folder like OUT.ABACUS
                natom = read_abacus_out.read_natom_fromlog(os.path.join(root, file))
                eks = read_abacus_out.read_e_fromlog(os.path.join(root, file))
                pressure = read_abacus_out.read_pressure_fromlog(os.path.join(root, file))
                bs = read_abacus_out.read_istate(os.path.join(root, "istate.info"))
                # continue if there is None among eks, pressure and bs
                parent = os.path.dirname(root)
                if None in [eks, pressure, bs]:
                    print(f"""WARNING: Present APNS job is broken: {parent}""")
                    continue
                bs = bs[0]
                atom_species, cellgen = read_apnsjob_desc(os.path.join(parent, "description.json"))
                system = os.path.basename(cellgen["config"])
                abacus_input = read_abacus_out.read_keyvals_frominput(os.path.join(parent, "INPUT"))
                ecutwfc = abacus_input["ecutwfc"]
                pps = [a["pp"] for a in atom_species]
                ppids = [": ".join(convert_fpp_to_ppid(pp)) for pp in pps]
                zvals = [float(ppparse.z_valence(os.path.join(parent, pp))) for pp in pps]
                s = "\n".join(ppids)
                print(f"""In folder {parent}
Structure tested: {system}
Number of atoms: {natom}
ecutwfc: {ecutwfc}
Final Kohn-Sham energy: {eks}
Pressure: {pressure}
Pseudopotentials are used:\n{s}
""")
                data = {"ecutwfc": ecutwfc, "eks": eks, "pressure": pressure, "istate": bs, "natom": natom, "z_valence": zvals}
                # band structure is not easy to print, therefore omitted
                idx = -1 if result.get(system) is None \
                    or result[system].get("ppcases") is None \
                        or result[system]["ppcases"].count(pps) == 0 \
                    else result[system]["ppcases"].index(pps)
                if idx == -1:
                    result.setdefault(system, {"ppcases": [], "pptests": []}).setdefault("ppcases", []).append(pps)
                    result[system]["pptests"].append([data])
                else:
                    result[system]["pptests"][idx].append(data)
                #result[(system, "|".join(pps), ecutwfc)] = (natom, zvals, eks, pressure, bs)
    return result

def repair_apnsjob(folder: str):
    """this function is written because one old abacustest version will drop description.json auto-generated by
    apns workflow. Now this has been fixed, so do not call this function in any cases.
    Call example:
    ```python
    repair_apnsjob("12310698")
    ```
    """
    import os
    import re
    import shutil
    for root, _, files in os.walk(folder):
        for file in files:
            if re.match(r"(running_)(\w+)(\.log)", file): # reach the inner most folder like OUT.ABACUS
                parent = os.path.dirname(root)
                f = os.path.basename(os.path.dirname(os.path.dirname(root)))
                shutil.copy2(os.path.join(f"../Yb_ecutwfc_test/{f}/description.json"), os.path.join(parent, "description.json"))

if __name__ == "__main__":
    from apns.analysis.apns2_ecut_utils import \
        update_ecutwfc, build_sptc_from_nested, plot_log, plot_stack
    import json, os
    from apns.analysis.apns2_utils import stru_rev_map
    sysrevmap_ = stru_rev_map("./apns_cache/structures.json", True)
    database = "normconserving-all-ecutwfctest"
    jobpath = "/root/documents/simulation/abacus/ultrasoft-test-4th-period"

    ###########################
    # update ecutwfc database #
    ###########################

    # collected = collect(jobpath)
    # system_and_stpcs = build_sptc_from_nested(collected)
    # result = []
    # for s, stpcs in system_and_stpcs.items():
    #     for stpc in stpcs:
    #         pp, data = stpc()
    #         temp = {"name": sysrevmap_[s], "fcif": s, "pp": pp}
    #         temp.update(data)
    #         result.append(temp)
    #         ecut_conv = stpc.ecuts[stpc.iconv]
    #         pp = stpc.pp(as_list=True)
    #         assert len(pp) == 1, "The pseudopotential should be unique for each test case"
    #         update_ecutwfc(pp[0], ecut_conv)

    # with open(os.path.basename(jobpath)+".json", "w") as f:
    #     json.dump(result, f)

    fdb = os.path.join("/root/abacus-develop/apns_toupdate", database + ".json")
    with open(fdb, "r") as f:
        result = json.load(f)

    ############################
    # plot ecutwfc convergence #
    ############################
    result = [r for r in result if r["name"] == "Si" and not r["pp"].endswith("(fr)")]
    result = sorted(result, key=lambda x: x["pp"])[:8]
    # flogs = plot_log(result)
    fstacks = plot_stack(result)
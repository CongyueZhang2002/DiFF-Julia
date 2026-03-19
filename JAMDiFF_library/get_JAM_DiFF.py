from pathlib import Path

import numpy as np

from tools.tools import load, lprint
from tools.config import conf, load_config
from fitlib.resman import RESMAN
from analysis.corelib import core

import warnings

warnings.filterwarnings("ignore")

def _resolve_wdir(wdir):
    wdir_path = Path(wdir)
    if wdir_path.is_absolute():
        return wdir_path
    if wdir_path.exists():
        return wdir_path.resolve()

    module_root = Path(__file__).resolve().parent
    candidate = module_root / wdir_path
    if candidate.exists():
        return candidate.resolve()

    raise FileNotFoundError(f"Could not resolve wdir: {wdir}")


def _mu_key(mu):
    return np.float16(mu)

def get_JAM_DiFF(z_array, M_array, Q2, kind, wdir):
    if kind not in {"H1a", "D1"} and wdir in {"H1a", "D1"}:
        kind, wdir = wdir, kind

    if kind == "H1a":
        conf_key = "tdiffpippim"
        method_name = "get_H"
        flavors = ["u"]
    elif kind == "D1":
        conf_key = "diffpippim"
        method_name = "get_D"
        flavors = ["u", "s", "c", "b", "g"]
    else:
        raise ValueError("kind must be 'H1a' or 'D1'")

    wdir = str(_resolve_wdir(wdir))

    print(f"\ngenerating {kind} at Q2 = {Q2} from {wdir}")
    load_config(f"{wdir}/input.py")
    istep = core.get_istep()
    replicas = core.get_replicas(wdir)
    core.mod_conf(istep, replicas[0])

    if Q2 is None:
        Q2 = conf["Q20"]

    conf["SofferBound"] = False
    resman = RESMAN(nworkers=1, parallel=False, datasets=False, load_lhapdf=False)
    parman = resman.parman
    parman.order = replicas[0]["order"][istep]

    jar = load(f"{wdir}/data/jar-{istep}.dat")
    replicas = jar["replicas"]

    diff_obj = conf[conf_key]
    Zgrid, Mgrid = np.meshgrid(z_array, M_array)
    Zgrid = Zgrid.flatten()
    Mgrid = Mgrid.flatten()
    L = len(Zgrid)

    out = {}
    for i, par in enumerate(replicas, start=1):
        lprint(f"{i}/{len(replicas)}")
        parman.set_new_params(par, initial=True)
        for flav in flavors:
            out.setdefault(flav, [])
            func = getattr(diff_obj, method_name)(
                Zgrid, Mgrid, Q2 * np.ones(L), flav, evolve=True
            )
            out[flav].append(func)

    return {"z_array": Zgrid, "Mh_array": Mgrid, "Q2": Q2, kind: out}

def get_JAM_DiFF_dict(z_array, Mh_array, mu_array, kinds, wdir):

    dict_raw_DiFF = {}
    keys = []
    flavors_by_kind = {}

    for kind in kinds:
        for mu in mu_array:
            obj = get_JAM_DiFF(z_array, Mh_array, mu**2, kind, wdir)
            mu_key = _mu_key(mu)

            value_dict = obj[kind]
            flavors = list(value_dict.keys())
            flavors_by_kind[kind] = flavors
            n_replicas = len(next(iter(value_dict.values())))

            for replica_id in range(n_replicas):
                key = (kind, mu_key, replica_id)
                keys.append(key)

                dict_raw_DiFF[key] = np.column_stack(
                    [value_dict[flav][replica_id] for flav in flavors]
                )

    dict_raw_DiFF["z_array"] = obj["z_array"]
    dict_raw_DiFF["Mh_array"] = obj["Mh_array"]
    dict_raw_DiFF["flavors_by_kind"] = flavors_by_kind
    if len(flavors_by_kind) == 1:
        dict_raw_DiFF["flavors"] = next(iter(flavors_by_kind.values()))

    return dict_raw_DiFF

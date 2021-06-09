from six import iteritems
from equilibrator_api import ComponentContribution
import numpy as np


api = ComponentContribution()


def fetch_metabolites(metabolites):
    """[summary]

    Parameters
    ----------
    metabolites : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    met_equilibrator_dict = {}
    for met, met_id in iteritems(metabolites):
        if met_id == "NA":
            eq_accession = None
        else:
            try:
                eq_accession = api.get_compound(met_id)
            except:
                eq_accession = None
        met_equilibrator_dict[met] = eq_accession

    return met_equilibrator_dict


def get_compound_vector(metabolite, met_id):
    equilibrator_access = fetch_metabolites({metabolite: met_id})
    equilibrator_instance = equilibrator_access[metabolite]
    if equilibrator_instance:
        try:
            rc_index = rc_compound_ids.index(equilibrator_instance.id)
            comp_vector = np.zeros(Nc + Ng, dtype=float)
            comp_vector[rc_index] = 1
            return comp_vector[np.newaxis, :]
        except ValueError:
            if equilibrator_instance.group_vector:
                comp_vector = np.hstack(
                    [
                        np.zeros(Nc, dtype=float),
                        equilibrator_instance.group_vector,
                    ]
                )
                return comp_vector[np.newaxis, :]
            else:
                comp_vector = np.zeros(Nc + Ng, dtype=float)
                return comp_vector[np.newaxis, :]
    else:
        comp_vector = np.zeros(Nc + Ng, dtype=float)
        return comp_vector[np.newaxis, :]

    def calculate_delG_f(metabolite, met_id, pH, I, temperature):
        """Calculates the standard transformed Gibbs formation energy of compound using component contribution method. pH, Ionic strength values are taken from model's compartment_info attribute

        Returns:
            float -- Transformed Gibbs energy of formation adjusted to pH, ionic strength of metabolite
        """
        compound_vector = get_compound_vector(metabolite, met_id)
        std_dG_f = compound_vector @ mu
        if self.compound_vector.any():
            transform = self.equilibrator_accession.transform(
                p_h=Q_(self.model.compartment_info["pH"][self.compartment]),
                ionic_strength=Q_(
                    str(self.model.compartment_info["I"][self.compartment]) + " M"
                ),
                temperature=Q_(str(default_T) + " K"),
            )
            return std_dG_f[0] + transform.to_base_units().magnitude * 1e-3
        else:
            return std_dG_f[0]
import os
import cobra
from cobra import Reaction, Metabolite

def get_or_create_met(model, met_id, name, compartment='c'):
    """Safely retrieves an existing metabolite or creates a new one."""
    try:
        return model.metabolites.get_by_id(met_id)
    except KeyError:
        new_met = Metabolite(met_id, name=name, compartment=compartment)
        model.add_metabolites([new_met])
        return new_met

def build_1_4_bdo_model():
    print("Building engineered iML1515 for 1,4-BDO (Paper 1)...")
    # Load base E. coli model directly from BiGG
    model = cobra.io.load_model("iML1515")
    
    # Existing cofactors and precursors in iML1515
    succoa = model.metabolites.get_by_id("succoa_c")
    sucsa = model.metabolites.get_by_id("sucsal_c") # Succinate semialdehyde
    coa = model.metabolites.get_by_id("coa_c")
    nadh = model.metabolites.get_by_id("nadh_c")
    nad = model.metabolites.get_by_id("nad_c")
    h = model.metabolites.get_by_id("h_c")
    
    # New engineered metabolites
    m_4hb = get_or_create_met(model, "M_4hb_c", "4-hydroxybutyrate")
    m_4hbcoa = get_or_create_met(model, "M_4hbcoa_c", "4-hydroxybutyryl-CoA")
    m_4hbald = get_or_create_met(model, "M_4hbald_c", "4-hydroxybutyraldehyde")
    m_14bdo = get_or_create_met(model, "M_14bdo_c", "1,4-butanediol")

    # Reaction 1: Succinyl-CoA reductase (SucD)
    # succoa + nadh + h -> sucsa + coa + nad
    rxn_sucD = Reaction("SUCD_eng")
    rxn_sucD.name = "CoA-dependent succinate semialdehyde dehydrogenase"
    rxn_sucD.add_metabolites({succoa: -1.0, nadh: -1.0, h: -1.0, sucsa: 1.0, coa: 1.0, nad: 1.0})
    
    # Reaction 2: 4-hydroxybutyrate dehydrogenase (4HBd)
    # sucsa + nadh + h -> 4hb + nad
    rxn_4hbd = Reaction("4HBD_eng")
    rxn_4hbd.add_metabolites({sucsa: -1.0, nadh: -1.0, h: -1.0, m_4hb: 1.0, nad: 1.0})

    # Reaction 3: 4-hydroxybutyryl-CoA transferase (Cat2)
    # 4hb + succoa -> 4hbcoa + succinate (Simplified to use generic CoA donor/acceptor for graph traversal)
    rxn_cat2 = Reaction("CAT2_eng")
    rxn_cat2.add_metabolites({m_4hb: -1.0, coa: -1.0, m_4hbcoa: 1.0})

    # Reaction 4: 4-hydroxybutyryl-CoA reductase (Ald)
    # 4hbcoa + nadh + h -> 4hbald + coa + nad
    rxn_ald = Reaction("ALD_eng")
    rxn_ald.add_metabolites({m_4hbcoa: -1.0, nadh: -1.0, h: -1.0, m_4hbald: 1.0, coa: 1.0, nad: 1.0})

    # Reaction 5: Alcohol dehydrogenase (Adh)
    # 4hbald + nadh + h -> 14bdo + nad
    rxn_adh = Reaction("ADH_eng")
    rxn_adh.add_metabolites({m_4hbald: -1.0, nadh: -1.0, h: -1.0, m_14bdo: 1.0, nad: 1.0})

    model.add_reactions([rxn_sucD, rxn_4hbd, rxn_cat2, rxn_ald, rxn_adh])
    
    # Save the model
    cobra.io.write_sbml_model(model, "models/e_coli_bdo_engineered.xml")
    print(" -> Saved as models/e_coli_bdo_engineered.xml\n")


def build_1_3_pdo_model():
    print("Building engineered iML1515 for 1,3-PDO (Paper 2)...")
    model = cobra.io.load_model("iML1515")
    
    # Existing cofactors and precursors
    glyc = model.metabolites.get_by_id("glyc_c") # Glycerol
    h2o = model.metabolites.get_by_id("h2o_c")
    nadh = model.metabolites.get_by_id("nadh_c")
    nad = model.metabolites.get_by_id("nad_c")
    h = model.metabolites.get_by_id("h_c")

    # New engineered metabolites
    m_3hpa = get_or_create_met(model, "M_3hpa_c", "3-hydroxypropionaldehyde")
    m_13pdo = get_or_create_met(model, "M_13pdo_c", "1,3-propanediol")

    # Reaction 1: Glycerol dehydratase (dhaB)
    # glycerol -> 3-hpa + h2o
    rxn_dhaB = Reaction("DHAB_eng")
    rxn_dhaB.add_metabolites({glyc: -1.0, m_3hpa: 1.0, h2o: 1.0})

    # Reaction 2: 1,3-propanediol oxidoreductase (dhaT)
    # 3-hpa + nadh + h -> 13pdo + nad
    rxn_dhaT = Reaction("DHAT_eng")
    rxn_dhaT.add_metabolites({m_3hpa: -1.0, nadh: -1.0, h: -1.0, m_13pdo: 1.0, nad: 1.0})

    model.add_reactions([rxn_dhaB, rxn_dhaT])
    
    cobra.io.write_sbml_model(model, "models/e_coli_pdo_engineered.xml")
    print(" -> Saved as models/e_coli_pdo_engineered.xml\n")


def build_artemisinin_model():
    print("Building engineered iMM904 for Artemisinin (Paper 3)...")
    model = cobra.io.load_model("iMM904")
    
    # Existing cofactors and precursors in yeast
    # Note: FPP is known as Farnesyl diphosphate (frdp_c) in iMM904
    fpp = model.metabolites.get_by_id("frdp_c") 
    nadph = model.metabolites.get_by_id("nadph_c")
    nadp = model.metabolites.get_by_id("nadp_c")
    nad = model.metabolites.get_by_id("nad_c")
    nadh = model.metabolites.get_by_id("nadh_c")
    h = model.metabolites.get_by_id("h_c")
    o2 = model.metabolites.get_by_id("o2_c")
    h2o = model.metabolites.get_by_id("h2o_c")

    # New engineered plant metabolites
    m_amph = get_or_create_met(model, "M_amph_c", "Amorphadiene")
    m_artalc = get_or_create_met(model, "M_artalc_c", "Artemisinic alcohol")
    m_artald = get_or_create_met(model, "M_artald_c", "Artemisinic aldehyde")
    m_artacid = get_or_create_met(model, "M_artemisinic_acid_c", "Artemisinic acid")

    # Reaction 1: Amorphadiene synthase (ADS)
    # fpp -> amorphadiene + ppi (ignoring ppi for simplicity in graph traversal)
    rxn_ads = Reaction("ADS_eng")
    rxn_ads.add_metabolites({fpp: -1.0, m_amph: 1.0})

    # Reaction 2: Cytochrome P450 (CYP71AV1) - Amorphadiene to artemisinic alcohol
    # amph + nadph + o2 + h -> artalc + nadp + h2o
    rxn_cyp1 = Reaction("CYP71AV1_eng")
    rxn_cyp1.add_metabolites({m_amph: -1.0, nadph: -1.0, o2: -1.0, h: -1.0, m_artalc: 1.0, nadp: 1.0, h2o: 1.0})

    # Reaction 3: Artemisinic alcohol dehydrogenase (ADH1)
    # artalc + nad -> artald + nadh + h
    rxn_adh1 = Reaction("ADH1_eng")
    rxn_adh1.add_metabolites({m_artalc: -1.0, nad: -1.0, m_artald: 1.0, nadh: 1.0, h: 1.0})

    # Reaction 4: Artemisinic aldehyde dehydrogenase (ALDH1)
    # artald + nad + h2o -> artacid + nadh + h
    rxn_aldh1 = Reaction("ALDH1_eng")
    rxn_aldh1.add_metabolites({m_artald: -1.0, nad: -1.0, h2o: -1.0, m_artacid: 1.0, nadh: 1.0, h: 1.0})

    model.add_reactions([rxn_ads, rxn_cyp1, rxn_adh1, rxn_aldh1])
    
    cobra.io.write_sbml_model(model, "models/s_cerevisiae_artemisinic_engineered.xml")
    print(" -> Saved as models/s_cerevisiae_artemisinic_engineered.xml\n")


if __name__ == "__main__":
    os.makedirs("models", exist_ok=True)
    build_1_4_bdo_model()
    build_1_3_pdo_model()
    build_artemisinin_model()
    print("All models successfully built! You can now run 'validate_papers.py'.")
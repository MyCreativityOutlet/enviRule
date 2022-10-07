package eawag.envirule;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;
import uk.ac.ebi.reactionblast.mechanism.MappingSolution;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class rule_generator {

    private boolean generalizeIgnoreHydrogen;
    private boolean includeFunctionalGroups;
    private String file;
    private int radius;

    private final SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
    private final SmilesGenerator sg = new SmilesGenerator(SmiFlavor.AtomAtomMap);
    private final Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.edgeShort()));

    public rule_generator(boolean generalizeIgnoreHydrogen, boolean includeFunctionalGroups, String file, int radius){
        this.generalizeIgnoreHydrogen = generalizeIgnoreHydrogen;
        this.includeFunctionalGroups = includeFunctionalGroups;
        this.file = file;
        this.radius = radius;
    }

    public Set<String> generate() throws  Exception {



        List<String> reactions = parseReactionFile(file);
        List<String> baseRules = new ArrayList<>();
        // Generalized form of rules
        List<String> newRules = new ArrayList<>();

        // Mapped reaction
        List<IReaction> atomMappingReactions = new ArrayList<>();

        // Final results
        Set<String> final_rules = new HashSet<>();

        // To cut off COA later
        int coa_count = 0;

        // Get rules for each reaction
        for(int i=0; i<reactions.size(); i++){

            // Unmapped form of each reaction
            String unmapped = reactions.get(i);
            String[] components = unmapped.split(">>"); // components[0]: reactants; components[1]: products
            String reactants = components[0];
            String products = components[1];

            // Remove catalysts from reaction
            // Note: now it can only remove catalysts if they show in both reactants and products.
            Set<String> reactants_set = new HashSet<>();
            Set<String> products_set = new HashSet<>();

            if(reactants.contains(".")){
                reactants_set.addAll(List.of(reactants.split("\\.")));
            }
            if(products.contains(".")){
                products_set.addAll(List.of(products.split("\\.")));
            }

            unmapped = "";

            if(reactants_set.size() != 0){
                for(String reactant: reactants.split("\\.")){
                    if(products_set.contains(reactant)) continue; // If one compound shows in both reactants and products
                    unmapped += reactant + ".";
                }
            }else{
                // It means there is only one reactant
                unmapped = reactants;
            }

            // If all compounds in reactants can also be found in products, then this reaction is problematic and should be skipped
            if(unmapped.length() == 0) continue;

            // Remove the redundant dot for reactants part
            if(unmapped.charAt(unmapped.length()-1) == '.'){
                unmapped = unmapped.substring(0, unmapped.length()-1);
            }

            // Append reaction sign
            unmapped += ">>";

            // Remove catalysts for products
            String unmapped_product = "";

            if(products_set.size() != 0){
                for(String product:products.split("\\.")){
                    if(reactants_set.contains(product)) continue;
                    unmapped_product += product + ".";
                }

                if(unmapped_product.charAt(unmapped_product.length()-1) == '.'){
                    unmapped_product = unmapped_product.substring(0, unmapped_product.length()-1);
                }
            }else{
                // If there is no dot in products, then there is only one single product
                unmapped_product += products;
            }

            // If everything in products can be found in reactants, then this reaction is problematic and should be skipped
            if(unmapped_product.length() == 0) continue;

            unmapped += unmapped_product;

            if(checkCoA(unmapped)) {
                coa_count += 1;
            }

            IReaction Reaction = sp.parseReactionSmiles(unmapped);
            IReaction performAtomAtomMapping = performAtomAtomMapping(Reaction, null);
            atomMappingReactions.add(performAtomAtomMapping);

            String base_rule = getNewRule(performAtomAtomMapping, sg, 0, false);









        }


        return null;

    }

    /**
     *
     * @param performAtomAtomMapping: atom-mapped reaction
     * @param sg
     * @param radius: radius to expand reaction center
     * @param unmapUnrelated: whether to set the mapping number of unmapped atoms to null
     * @return
     * @throws Exception
     */
    private String getNewRule(IReaction performAtomAtomMapping, SmilesGenerator sg, int radius, Boolean unmapUnrelated) throws Exception {

        // Get changed atoms
        Set<Integer> changed_atom_tags = getChangedAtoms(performAtomAtomMapping, radius);

        // Extract reaction center from reactants and products
        IAtomContainer reactant_fragments = getFragmentsMol(performAtomAtomMapping.getReactants(), changed_atom_tags);
        IAtomContainer product_fragments = getFragmentsMol(performAtomAtomMapping.getProducts(), changed_atom_tags);

        return cleanSMIRKS(reactant_fragments, product_fragments, sg, changed_atom_tags, performAtomAtomMapping, unmapUnrelated);
    }

    private Set<Integer> getChangedAtoms(IReaction performAtomAtomMapping, int radius) throws Exception{
        // Get changed atoms
        BondChangeCalculator bcc = new BondChangeCalculator(performAtomAtomMapping);
        bcc.computeBondChanges(true, false);
        Set<Integer> changed_atom_tags = new HashSet<Integer>();

        for(IAtom atom: bcc.getReactionCenterSet()){
            changed_atom_tags.add(atom.getMapIdx());
        }

        // expand the neighbors in reactants
        for(int i=0; i<radius; i++){
            expandReactantNeighbors(performAtomAtomMapping.getReactants(), changed_atom_tags);
        }

        // include functional groups
        if(includeFunctionalGroups){
            if(radius != 0) addFunctionalGroups(performAtomAtomMapping.getReactants(), changed_atom_tags);
        }

        // include all the added groups in products if there are any
        newAddedGroupsInProducts(performAtomAtomMapping, changed_atom_tags);

        return changed_atom_tags;
    }

    /**
     * For all the atoms only in products, they should be added to the reaction center.
     * @param performAtomAtomMapping
     * @param changed_atom_tags
     */
    private void newAddedGroupsInProducts(IReaction performAtomAtomMapping, Set<Integer> changed_atom_tags){

        IAtomContainerSet reactants = performAtomAtomMapping.getReactants();
        IAtomContainerSet products = performAtomAtomMapping.getProducts();

        Set<Integer> reactants_atom_tags = new HashSet<>();

        for(int i=0; i<reactants.getAtomContainerCount(); i++){
            IAtomContainer mol = reactants.getAtomContainer(i);
            for(IAtom atom: mol.atoms()){
                reactants_atom_tags.add(atom.getMapIdx());
            }
        }

        for(int i=0; i<products.getAtomContainerCount(); i++){
            IAtomContainer mol = products.getAtomContainer(i);
            for(IAtom atom: mol.atoms()){
                if(!reactants_atom_tags.contains(atom.getMapIdx())){
                    changed_atom_tags.add(atom.getMapIdx());
                }
            }
        }
    }

    private void addFunctionalGroups(IAtomContainerSet reactants, Set<Integer> changed_atom_tags){
        // The list of functional groups are from
        // "Prediction of Organic Reaction Outcomes Using Machine Learning"
        // https://github.com/connorcoley/ochem_predict_nn/blob/master/data/generate_reaction_templates.py

        Set<Integer> addedFunctionalGroups = new HashSet<>();

        String[] group_templates = {
                "C(=O)Cl", //acid chloride
                "C(=O)O", // carboxylic acid
                "C(=O)[O-]", // carboxylic acid
                "[S](=O)(=O)(Cl)", //sulfonyl chloride
                "[B](O)(O)", //boronic acid
                "[N](=!@C=!@O)", //isocyanate
                "[N]=[N]=[N]", //azide
                "O=C1N(Br)C(=O)CC1", //NBS brominating agent
                "C=O", //carbonyl
                "ClS(Cl)=O", //thionyl chloride
                "[Mg][Br,Cl]", //grinard(non - disassociated)
                "[#6]S(=O)(=O)[O]", //RSO3 leaving group
                "[O]S(=O)(=O)[O]", //SO4 group
                "[N-]=[N+]=[C]", //diazo - alkyl
        };

        for(int i=0; i<reactants.getAtomContainerCount(); i++){
            IAtomContainer mol = reactants.getAtomContainer(i);
            for(String smarts: group_templates){
                Pattern pattern = SmartsPattern.create(smarts);
                Mappings matches = pattern.matchAll(mol);

                if(matches.count() == 0) continue;
                // matches are int[][] of atom index (not the MapId of atom).
                for(int[] match: matches){
                    Set<Integer> group_match = new HashSet<>();
                    for(int m: match) group_match.add(mol.getAtom(m).getMapIdx());

                    for(IAtom atom: mol.atoms()){
                        // Make sure the atom in reactant is in reaction center and also in a functional group.
                        if(changed_atom_tags.contains(atom.getMapIdx()) && group_match.contains(atom.getMapIdx())){
                            addedFunctionalGroups.addAll(group_match);
                        }
                    }
                }
            }
        }
        changed_atom_tags.addAll(addedFunctionalGroups);
    }

    private void expandReactantNeighbors(IAtomContainerSet compounds, Set<Integer> changed_atom_tags){

        Set<Integer> add_atom_tags = new HashSet<Integer>();

        for(int i=0; i<compounds.getAtomContainerCount(); i++){
            IAtomContainer mol = compounds.getAtomContainer(i);
            for(IBond bond: mol.bonds()){

                // For a bond, if one end atom (A) is included in reaction center, while the other end atom (B) is not. Then B should be added into reaction center.
                int count = 0;
                if(changed_atom_tags.contains(bond.getEnd().getMapIdx())) count ++;
                if(changed_atom_tags.contains(bond.getBegin().getMapIdx())) count ++;
                if(count == 1){
                    add_atom_tags.add(bond.getBegin().getMapIdx());
                    add_atom_tags.add(bond.getEnd().getMapIdx());
                }
            }
        }

        changed_atom_tags.addAll(add_atom_tags);
    }


    private IReaction performAtomAtomMapping(IReaction cdkReaction, String reactionName) throws InvalidSmilesException, AssertionError, Exception {
        cdkReaction.setID(reactionName);
        /*
         RMT for the reaction mapping
         */
        boolean forceMapping = true;//Overrides any mapping present int the reaction
        boolean generate2D = true;//2D perception of the stereo centers
        boolean generate3D = false;//2D perception of the stereo centers
        boolean standardizeReaction = true; //Standardize the reaction
        ReactionMechanismTool rmt = new ReactionMechanismTool(cdkReaction, forceMapping, generate2D, generate3D, standardizeReaction);
        MappingSolution s = rmt.getSelectedSolution();//Fetch the AAM Solution
        IReaction reaction = s.getReaction();//Fetch Mapped Reaction
        return reaction;
    }

    private boolean checkCoA(String base_rule) throws CDKException {

        Pattern ptrn1 = SmartsPattern.create("CC(C)(COP(O)(=O)OP(O)(=O)OCC1OC(C(O)C1OP(O)(O)=O)n1cnc2c(N)ncnc12)C(O)C(=O)NCCC(=O)NCCS");
        String reactants = base_rule.split(">>")[0];
        String products = base_rule.split(">>")[1];

        IAtomContainer reactants_mol = sp.parseSmiles(reactants);
        IAtomContainer products_mol = sp.parseSmiles(products);

        aromaticity.apply(reactants_mol);
        aromaticity.apply(products_mol);

        Mappings reactants_res = ptrn1.matchAll(reactants_mol);
        Mappings products_res = ptrn1.matchAll(products_mol);

        if (reactants_res.count() == 0 && products_res.count() != 0){
            return true;
        }

        return false;
    }

    private List<String> parseReactionFile(String file) throws IOException {
        List<String> reactions = new ArrayList<>();

        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        while ((line = br.readLine()) != null) {
            reactions.add(line);
        }
        br.close();
        return reactions;
    }


}

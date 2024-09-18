package eawag.envirule.modules;

import ambit2.core.data.MoleculeTools;
import ambit2.smarts.SMIRKSManager;
import ambit2.smarts.SMIRKSReaction;
import ambit2.smarts.SmartsConst;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.parser.SMARTSParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class reactor {

    private static final LinkedHashMap<String, String> flattenRegs = new LinkedHashMap<>();
    static {
        flattenRegs.put("@", "");
    }
    private static final LinkedHashMap<String, String> uncistransRegs = new LinkedHashMap<>();
    static {
        uncistransRegs.put("/", "");
        uncistransRegs.put("\\\\", "");
    }
    private static final Set<String> basicRuleSet = new HashSet<>();
    static {
        basicRuleSet.add("[H][N+:1]([H])([H])[#6:2]>>[H][#7:1]([H])-[#6:2]");
        basicRuleSet.add("[H][#8:1][C:2]#[N:3]>>[#8-:1][C:2]#[N:3]");
        basicRuleSet.add("[H][#8:1]-[#6:2]=[O:3]>>[#8-:1]-[#6:2]=[O:3]");
        basicRuleSet.add("[H][#8:1]-[#7+:2](-[*:3])=[O:4]>>[#8-:1]-[#7+:2](-[*:3])=[O:4]");
        basicRuleSet.add("[#6;A:1][#6:2](-[#8-:3])=[#6;A:4]>>[#6:1]-[#6:2](-[#8:3][H])=[#6;A:4]");
        basicRuleSet.add("[H][#8:1]-[$([#15]);!$(P([O-])):2]>>[#8-:1]-[#15:2]");
        basicRuleSet.add("[H][#8:1]-[c:2]1[c:3][c:4][c:5]([c:6][c:7]1-[#7+:8](-[#8-:9])=[O:10])-[#7+:11](-[#8-:12])=[O:13]>>[#8-:1]-[c:2]1[c:3][c:4][c:5]([c:6][c:7]1-[#7+:8](-[#8-:9])=[O:10])-[#7+:11](-[#8-:12])=[O:13]");
        basicRuleSet.add("[H][#8:1][S:2]([#8:3][H])(=[O:4])=[O:5]>>[#8-:1][S:2]([#8-:3])(=[O:4])=[O:5]");
        basicRuleSet.add("[#6:1]-[#8:2][S:3]([#8:4][H])(=[O:5])=[O:6]>>[#6:1]-[#8:2][S:3]([#8-:4])(=[O:5])=[O:6]");
        basicRuleSet.add("[H][#8:3][S:2]([#6:1])(=[O:4])=[O:5]>>[#6:1][S:2]([#8-:3])(=[O:4])=[O:5]");
        basicRuleSet.add("[H][#8:1][S:2]([*:3])=[O:4]>>[#8-:1][S:2]([*:3])=[O:4]");
        basicRuleSet.add("[N+:1]([H])>>[N:1]");
        basicRuleSet.add("[H][#8:1]-[c:2]>>[#8-:1]-[c:2]");
    }

    private static final Set<String> enhancedRuleSet = new HashSet<>();
    static {
        enhancedRuleSet.addAll(basicRuleSet);
        enhancedRuleSet.add("[H][#8:1]-[#15:2]>>[#8-:1]-[#15:2]");
    }

    private static final Set<String> exoticRuleSet = new HashSet<>();
    static {
        exoticRuleSet.addAll(enhancedRuleSet);
        exoticRuleSet.add("[H][S:1]-[#15:2]=[$([#16]),$([#8]):3]>>[S-:1]-[#15:2]=[$([#16]),$([#8]):3]");
    }

    private static final Set<String> cutCoARuleSet = new HashSet<>();
    static {
        cutCoARuleSet.addAll(exoticRuleSet);
        cutCoARuleSet.add("CC(C)(COP(O)(=O)OP(O)(=O)OCC1OC(C(O)C1OP(O)(O)=O)n1cnc2c(N)ncnc12)C(O)C(=O)NCCC(=O)NCCS[$(*):1]>>[O-][$(*):1]");
    }

    private static final Set<String> enolKetoRuleSet = new HashSet<>();
    static {
        enolKetoRuleSet.addAll(cutCoARuleSet);
        enolKetoRuleSet.add("[H][#8:2]-[#6:3]=[#6:1]>>[#6:1]-[#6:3]=[O:2]");
    }

    public static final Set<String> allStandardizeRules = new HashSet<>();
    static{
        allStandardizeRules.addAll(enolKetoRuleSet);
    }

    public static final Set<LinkedHashMap<String, String>> allStandardizeRegRules = new HashSet<>();
    static{
        allStandardizeRegRules.add(flattenRegs);
        allStandardizeRegRules.add(uncistransRegs);
    }

    public static String standardizeSmiles(String smiles, Set<String> standardizeRules, Set<LinkedHashMap<String, String>> standardizeRegRules) throws Exception {

        String normalizedSmiles = toSmiles(fromSmiles(smiles));

        for (LinkedHashMap<String, String> rule: standardizeRegRules){
            normalizedSmiles = RegxStandardizeSmiles(rule, normalizedSmiles);
        }

        for (String rule: standardizeRules){
            normalizedSmiles = RuleStandardizeSmiles(rule, normalizedSmiles);
        }

        return normalizedSmiles;
    }

    public static class CommandLineParser {
        List <String> args;
        HashMap<String, List<String>> map = new HashMap<>();
        Set<String> flags = new HashSet<>();

        CommandLineParser(String arguments[])
        {
            this.args = Arrays.asList(arguments);
            map();
        }

        // Return argument names
        public Set<String> getArgumentNames()
        {
            Set<String> argumentNames = new HashSet<>();
            argumentNames.addAll(flags);
            argumentNames.addAll(map.keySet());
            return argumentNames;
        }

        // Check if flag is given
        public boolean getFlag(String flagName)
        {
            return flags.contains(flagName);
        }

        // Return argument value for particular argument name
        public String[] getArgumentValue(String argumentName)
        {
            if(map.containsKey(argumentName))
                return map.get(argumentName).toArray(new String[0]);
            else
                return null;
        }

        // Map the flags and argument names with the values
        public void map()
        {
            for(String arg: args)
            {
                if(arg.startsWith("-"))
                {
                    if (args.indexOf(arg) == (args.size() - 1))
                    {
                        flags.add(arg.replace("-", ""));
                    }
                    else if (args.get(args.indexOf(arg)+1).startsWith("-"))
                    {
                        flags.add(arg.replace("-", ""));
                    }
                    else
                    {
                        //List of values (can be multiple)
                        List<String> argumentValues = new ArrayList<>();
                        int i = 1;
                        while(args.indexOf(arg)+i != args.size() && !args.get(args.indexOf(arg)+i).startsWith("-"))
                        {
                            argumentValues.add(args.get(args.indexOf(arg)+i));
                            i++;
                        }
                        map.put(arg.replace("-", ""), argumentValues);
                    }
                }
            }
        }
    }

    public static String RuleStandardizeSmiles(String rule, String smiles) throws Exception{

        List<IAtomContainer[]> products;
        products = apply(preprocessMolecule(smiles), false, false, rule);
        if (products.isEmpty()) return smiles;

        if(products.size() > 1 || products.get(0).length > 1){
            String msg = "ERROR in applyRepeatedly: >1 results for smiles '"+smiles+"' and rule '"+rule+":\n";
            for (IAtomContainer[] iAtomContainers : products) {
                for (IAtomContainer iAtomContainer : iAtomContainers) {
                    msg += toSmiles(iAtomContainer)+" ";
                }
                msg+="\n";
            }
            // System.err.println(msg);
        }

        return toSmiles(products.get(0)[0]);
    }

    public static String RegxStandardizeSmiles(LinkedHashMap<String, String> rule, String smiles) throws Exception {
        String modified = smiles;
        for (Map.Entry<String, String> substitute : rule.entrySet())
            modified  = modified.replaceAll(substitute.getKey(),substitute.getValue());

        return toSmiles(fromSmiles(modified));
    }

    public static void main(String[] args) throws Exception {

        CommandLineParser clp = new CommandLineParser(args);

        if(!clp.map.containsKey("rules") || clp.getArgumentValue("rules").length != 1){
            System.out.println("One and only one rule file please.");
            System.exit(1);
        }

        if(!clp.map.containsKey("reactions") || clp.getArgumentValue("reactions").length != 1){
            System.out.println("One and only one reaction file please.");
            System.exit(1);
        }

        boolean manual = clp.flags.contains("manual");

        String reactions_file = clp.getArgumentValue("reactions")[0];
        String rules_file = clp.getArgumentValue("rules")[0];

        List<String> smirks_list = parseRuleFile(rules_file);
        List<String> reactions = parseReactionFile(reactions_file);

        for (String reaction : reactions) {

            String reactant = reaction.split(">>")[0];

            IAtomContainer Reactant = preprocessMolecule(reactant);
            System.out.println("Reactant: ");
            System.out.println(toSmiles(Reactant));

            // Atoms will be perceived in the method "validateMolecule"
            if (!manual) {
                validateMolecule(Reactant);
            }

            List<IAtomContainer[]> ProductLists;
            try {
                ProductLists = apply_rules(smirks_list, Reactant);
            } catch (Exception e) {
                System.out.println("Exception on " + reactant);
                continue;
            }

            if (ProductLists.size() == 0) System.out.println("No products for: " + reactant);

            Set<String> ruleProductSmiles = new HashSet<>();

            for (IAtomContainer[] productList : ProductLists) {
                for (IAtomContainer ruleproduct : productList) {

                    String productSmiles = toSmiles(ruleproduct);
                    IAtomContainer mol = preprocessMolecule(productSmiles);

                    productSmiles = toSmiles(mol);

                    productSmiles = standardizeSmiles(productSmiles, allStandardizeRules, allStandardizeRegRules);
                    ruleProductSmiles.add(productSmiles);
                }
            }
            System.out.println("Products:");
            for (String product_smiles : ruleProductSmiles)
                System.out.println(product_smiles);
        }
    }

    public static List<String> parseRuleFile(String file) throws IOException{
        List<String> rules = new ArrayList<>();

        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        while ((line = br.readLine()) != null) {
            rules.add(line);
        }

        br.close();

        return rules;
    }

    private static List<IAtomContainer[]> apply_rules(List<String> rule_list, IAtomContainer preprocessedMolecule) throws Exception {

        List<IAtomContainer[]> list = new ArrayList<>();
        // init result list with mol
        list.add(new IAtomContainer[]{preprocessedMolecule});
        List<IAtomContainer[]> newList = new ArrayList<>();
        for (String rule : rule_list) {
            if(checkSmirks(rule)){
                for (IAtomContainer[] substrates : list) {
                    for (IAtomContainer substrate : substrates) {
                        List<IAtomContainer[]> products = apply(substrate, true, true, rule);
                        newList.addAll(products);
                    }
                }
            }
        }
        return newList;
    }

    public static boolean checkSmirks(final String newSmirks) throws Exception {
        SMIRKSReaction reaction = fromSmirks(newSmirks);

        if (reaction==null) {
            checkValidSmartsException(newSmirks);
            return false;
        }

        return true;
    }

    protected static void checkValidSmartsException(String smarts)
            throws Exception {
        try {

            QueryAtomContainer q = SMARTSParser
                    .parse(smarts, SilentChemObjectBuilder.getInstance());
            if (q == null) {
                throw new NullPointerException();
            }
        } catch (Throwable e) {
            throw new Exception("not a valid smarts string '"
                    + smarts + "', error: " + e.getMessage());
        }
    }

    public static List<IAtomContainer[]> apply(IAtomContainer preprocessedMolecule,
                                                boolean splitUnconnectedMolecules,
                                                boolean applyRuleToSinglePosition,
                                                String smirks)
            throws Exception {
        List<IAtomContainer[]> rProducts = new ArrayList<>();
        HashSet<String> rProductSmiles = new HashSet<>();
        try {

            IAtomContainerSet products = applySmirks(
                    smirks,
                    preprocessedMolecule,
                    applyRuleToSinglePosition);


            if (products != null) {

                for (IAtomContainer p : products.atomContainers()) {
                    if (splitUnconnectedMolecules) {
                        List<IAtomContainer> splitProducts = new
                                ArrayList<>();
                        for (IAtomContainer m : ConnectivityChecker
                                .partitionIntoMolecules(p).atomContainers()) {

                            String smiles = toSmiles(m);
                            if (!rProductSmiles.contains(smiles)) {
                                rProductSmiles.add(smiles);
                                splitProducts.add(m);
                            }
                        }
                        rProducts.add(splitProducts.toArray(new IAtomContainer[splitProducts.size()]));
                    } else {
                        String smiles = toSmiles(p);
                        if (!rProductSmiles.contains(smiles)) {
                            rProductSmiles.add(smiles);
                            rProducts.add(new IAtomContainer[]{p});
                        }
                    }
                }
            }
        } catch (Exception e) {
            throw new Exception("Error while applying Rule to compound '" + preprocessedMolecule
                    .toString()
                    + "'.", e);
        }
        return rProducts;
    }

    public static IAtomContainerSet applySmirks(final String smrk,
                                                final IAtomContainer
                                                        preprocessedMolecule,
                                                final boolean singlePos) throws
            Exception {
        SMIRKSManager smrkMan = new SMIRKSManager(SilentChemObjectBuilder
                .getInstance());
        smrkMan.setFlagSSMode(SmartsConst.SSM_MODE.SSM_NON_IDENTICAL_FIRST);
        smrkMan.setFlagProcessResultStructures(true);
        smrkMan.setFlagClearHybridizationBeforeResultProcess(true);
        smrkMan.setFlagClearImplicitHAtomsBeforeResultProcess(true);
        smrkMan.setFlagClearAromaticityBeforeResultProcess(true);
        smrkMan.setFlagAddImplicitHAtomsOnResultProcess(true);
        smrkMan.setFlagConvertAddedImplicitHToExplicitOnResultProcess(false);
        smrkMan.setFlagConvertExplicitHToImplicitOnResultProcess(true);
        smrkMan.getSmartsParser().mSupportDoubleBondAromaticityNotSpecified =
                false;
        smrkMan.setFlagApplyStereoTransformation(true);

        SMIRKSReaction reaction = smrkMan.parse(smrk);
        if (!smrkMan.getErrors().equals(""))
            throw new Exception("Invalid SMIRKS: " + smrkMan.getErrors());

        IAtomContainerSet set;
        if (singlePos)
            set = smrkMan.applyTransformationWithSingleCopyForEachPos(
                    preprocessedMolecule, null, reaction,
                    SmartsConst.SSM_MODE.SSM_ALL);
        else
            set = smrkMan.applyTransformationWithCombinedOverlappedPos(
                    preprocessedMolecule, null, reaction);

        if (set != null) {
            IAtomContainerSet postProcessedMols = new AtomContainerSet();
            for (IAtomContainer mol : set.atomContainers()) {
                mol = AtomContainerManipulator.suppressHydrogens(mol);
                postProcessedMols.addAtomContainer(mol);
            }
            set = postProcessedMols;
        }
        return set;
    }

    public static void validateMolecule(final IAtomContainer molecule)
            throws CDKException {
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
        for (int i = 0; i < molecule.getAtomCount(); i++) {
            if (null == molecule.getAtom(i).getValency()) {
                throw new CDKException("Valency is not available");
            }
        }

    }

    public static SMIRKSReaction fromSmirks(final String smirks)
            throws Exception {
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        final SMIRKSManager ms = new SMIRKSManager(builder);
        SMIRKSReaction react;
        try {
            react = ms.parse(smirks);
        } catch (Exception e) {
            throw new Exception("Error parsing SMIRKS '" + smirks +
                    "': "
                    + e.getMessage(), e);
        }
        if (!ms.getErrors().equals("")) {
            throw new Exception("Error parsing SMIRKS '" + smirks
                    + "': " + ms.getErrors());
        }
        return react;
    }

    public static String toSmiles(final IAtomContainer mol)
            throws Exception {
        try {
            return SmilesGenerator.absolute().create(AtomContainerManipulator
                    .copyAndSuppressedHydrogens(mol));
        } catch (final Exception e) {
            throw new Exception("Molecule: " + mol, e);
        }
    }

    public static IAtomContainer fromSmiles(final String smiles)
            throws Exception {

        // cdk1.5 start
        try {
            final SmilesParser sp = new SmilesParser(
                    SilentChemObjectBuilder.getInstance());

            return sp.parseSmiles(smiles);
        } catch (final Exception e) {
            throw new Exception("SMILES: " + smiles, e);
        }
        // cdk1.5 end

    }

    public static IAtomContainer preprocessMolecule(final String smiles)
            throws Exception {
        try {
            IAtomContainer target = fromSmiles(smiles);
            for (IAtom atom : target.atoms())
                if (atom.getFlag(CDKConstants.ISAROMATIC))
                    atom.setFlag(CDKConstants.ISAROMATIC, false);
            for (IBond bond : target.bonds())
                if (bond.getFlag(CDKConstants.ISAROMATIC))
                    bond.setFlag(CDKConstants.ISAROMATIC, false);

            MoleculeTools.convertImplicitToExplicitHydrogens(target);

            Aromaticity aromaticity = new Aromaticity(
                    ElectronDonation.daylight(),
                    Cycles.or(Cycles.all(), Cycles.edgeShort()));
            aromaticity.apply(target);
            return target;
        } catch (Exception e) {
            throw new Exception("Cannot preprocess compound", e);
        }
    }

    public static List<String> parseReactionFile(String file) throws IOException {
        List<String> reactions = new ArrayList<>();

        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        while ((line = br.readLine()) != null) {

            // Remove catalysts
            String[] s = line.split(">>");

            if(s[0].contains(".")){
                String[] reactants = s[0].split("\\.");
                String[] products = s[1].split("\\.");
                Set<String> reactants_set = new HashSet<>(Arrays.asList(reactants));
                Set<String> products_set = new HashSet<>(Arrays.asList(products));
                if(reactants_set.size() != products_set.size()) continue;
                line = "";
                for(String reactant:reactants){
                    if(products_set.contains(reactant)) continue;
                    line += reactant + ".";
                }
                if(line.length() == 0) continue;
                if(line.charAt(line.length()-1) == '.'){
                    line = line.substring(0, line.length()-1);
                }
                line += ">>";
                for(String product:products){
                    if(reactants_set.contains(product)) continue;
                    line += product + ".";
                }
                if(line.charAt(line.length()-1) == '.'){
                    line = line.substring(0, line.length()-1);
                }
            }

            reactions.add(line);
        }

        br.close();
        return reactions;
    }


}

package eawag.envirule;

import eawag.envirule.modules.rule_generator;
import eawag.envirule.modules.reactor;
import org.openscience.cdk.interfaces.IAtomContainer;
import py4j.GatewayServer;


import java.util.*;
import java.util.List;

import static eawag.envirule.modules.reactor.parseReactionFile;
import static eawag.envirule.modules.reactor.standardizeSmiles;
import static eawag.envirule.modules.rxn_clusterer.*;

public class enviRule_demo {
    public static Set<String> generateRules(String file, boolean generalizeIgnoreHydrogen, boolean includeFunctionalGroups, int radius) throws Exception {

        Map<String, Map<Integer, Double>> reactionFPs = new HashMap<String, Map<Integer, Double>>();
        List<String> reactions = parseReactions(file, smilesParser, reactionFPs);
        Map<String, Set<String>> clusters = getClusters(reactions, smilesParser, reactionFPs, "out_dir", false);
        System.out.println("Found " + String.valueOf(clusters.size()) + " reaction clusters");
        Set<String> rules = new HashSet<>();
        for (Map.Entry<String, Set<String>> entry: clusters.entrySet()) {
            List<String> cluster = new ArrayList<>(entry.getValue());
            if (cluster.size() >= 2) {
                rule_generator generator = new rule_generator(generalizeIgnoreHydrogen, includeFunctionalGroups, cluster, radius);
                rules.addAll(generator.generate());
            }
        }
        System.out.println("Finished generating rules, found " + String.valueOf(rules.size()));
        return rules;
    }

    public static String standardSmiles(String smiles) {
        try {
            return reactor.standardizeSmiles(smiles, reactor.allStandardizeRules, reactor.allStandardizeRegRules);
        }
        catch (Exception e) {
            return "";
        }
    }

    public static String runReaction(String smiles, String rule) {
        try {
            IAtomContainer reactant = reactor.preprocessMolecule(smiles);
            List<IAtomContainer[]> all_products = new ArrayList<>();
            if (reactor.checkSmirks(rule)) {
                List<IAtomContainer[]> products = reactor.apply(reactant, true, true, rule);
                all_products.addAll(products);
            }
            Set<String> ruleProductSmiles = new HashSet<>();

            for (IAtomContainer[] productList : all_products) {
                for (IAtomContainer ruleproduct : productList) {

                    String productSmiles = reactor.toSmiles(ruleproduct);
                    IAtomContainer mol = reactor.preprocessMolecule(productSmiles);

                    productSmiles = reactor.toSmiles(mol);

                    productSmiles = reactor.standardizeSmiles(productSmiles, reactor.allStandardizeRules, reactor.allStandardizeRegRules);
                    ruleProductSmiles.add(productSmiles);
                }
            }
            return String.join(".", ruleProductSmiles);
        } catch (Exception e) {
            return "";
        }
    }

    public enviRule_demo() {

    }

    public static void main(String[] args) throws Exception {
        boolean gateway = args.length == 0;
        if (gateway) {
            GatewayServer gatewayServer = new GatewayServer(new enviRule_demo());
            gatewayServer.start();
            System.out.println("Gateway Server Started");
        } else {
            boolean generalizeIgnoreHydrogen = true;
            boolean includeFunctionalGroups = true;
            String file = "example/temp_reactions.txt";
            List<String> reactions = parseReactionFile(file);
            System.out.println("Testing standardisation");
            String[] components = reactions.get(0).split(">>"); // components[0]: reactants; components[1]: products
            String standardised = "";
            for (int i = 0; i < 1; i++) {
                String reactants = "CCC(CCCC(C(C)C(=O)C1=CC2C(C=C(C)C3CC(CC32)OC(=C(C(C(C(=O)C)OC)OCC)O)O)C1CC(=O)[O-])OC4CCC(C(C)O4)N(C)C)[O-]";
                standardised = enviRule_demo.standardSmiles(reactants);
                System.out.println("Standardised " + reactants + " to " + standardised);
            }

            int radius = 1;
            System.out.println("Testing rule generation");
            Set<String> rules = enviRule_demo.generateRules(file, true, true, 3);

            System.out.println("Finished rule generation");
            System.out.println("Testing reaction running");
            for (String rule: rules) {
                enviRule_demo.runReaction(standardised, rule);
                break;
            }
            for (String rule: rules) {
                System.out.println(rule);
            }
            System.exit(0);
        }
    }

}

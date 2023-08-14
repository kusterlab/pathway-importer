"""Functions to download and convert database dumps from various Pathway databases.
The raw files are always of some XML format, the processed files are in a unified JSON format that can be parsed
by D3.js."""

import re
import requests
import json
import zipfile
import pandas as pd
from io import BytesIO
from collections import defaultdict
from typing import Union, TypedDict
from time import sleep
from pathlib import Path
from lxml import etree
from tqdm import tqdm
import constants


def get_kegg_pathway_list(organism_identifier: str) -> pd.DataFrame:
    print('Retrieving pathway identifiers...')
    pathway_list_response = requests.get(constants.KEGG_PATHWAY_LIST_ENDPOINT.format(organism_identifier))
    pathwayid_to_pathwayname_df = pd.read_csv(BytesIO(pathway_list_response.content), sep='\t',
                                              names=['pathway_id', 'name'])
    # We actually only need the pathway_ids
    pathway_ids = pathwayid_to_pathwayname_df.pathway_id
    return pathway_ids


def download_kegg_pathways(target_dir: Path, organism_identifier: str) -> None:
    target_dir.mkdir(exist_ok=True, parents=True)
    pathway_ids = get_kegg_pathway_list(organism_identifier)
    for pathway_id in tqdm(pathway_ids, 'Downloading KEGG Pathways'):
        raw_kgml_outputfile = Path(target_dir, '{}.kgml'.format(pathway_id))
        kgml_response = requests.get(constants.KEGG_PATHWAY_KGML_GENERIC_ENDPOINT.format(pathway_id))
        if len(kgml_response.content) == 0:
            print('\tKGML file not found: {}'.format(pathway_id))
            continue
        with open(raw_kgml_outputfile, 'wb') as o:
            o.write(kgml_response.content)


def download_wikipathways_pathways(target_dir: Path, organism_identifier: str) -> None:
    print('Downloading latest Wikipathways Dump...')
    target_dir.mkdir(exist_ok=True, parents=True)
    # Get date of current wikipathways dump
    wp_current_gpml_html = requests.get('https://wikipathways-data.wmcloud.org/current/gpml/')
    last_database_dump_date = re.search('wikipathways-\d+', wp_current_gpml_html.text).group().split('-')[-1]
    zip_response = requests.get(
        constants.WIKIPATHWAYS_ALL_PATHWAYS_ZIPPED_ENDPOINTS_MAP[organism_identifier].format(last_database_dump_date),
        stream=True)
    zip_file = zipfile.ZipFile(BytesIO(zip_response.content))
    print('Extracting Wikipathways Pathways...')
    zip_file.extractall(target_dir)
    print('Done.')


def map_gene_identifiers(from_database: str, to_database: str, gene_ids: set) -> defaultdict[list]:
    payload = {
        'from': from_database,
        'to': to_database,  # 'UniProtKB'
        'ids': ','.join(gene_ids)
    }
    request = requests.post(constants.UNIPROT_MAPPING_ENDPOINT, payload)
    jobid = json.loads(request.text).get('jobId')
    job_finished = False
    while not job_finished:
        result_response = json.loads(requests.get(constants.UNIPROT_RESULT_ENDPOINT + jobid).text)
        job_finished = (result_response.get('results') is not None)
        if not job_finished:
            print('Job pending: {}'.format(jobid))
        # Be nice to Uniprot and wait a second before asking again
        sleep(1)

    query_result = defaultdict(list)
    for query in result_response['results']:
        query_result[query['from']].append(query['to'])
    return query_result


def map_ec_to_uniprot(taxcode: int, ec_numbers: set) -> dict[list]:
    """For EC Identifiers, we cannot use the UNIPROT_MAPPING_ENDPOINT as in the function above
    so they are queried somewhat differently"""
    ec_to_uniprot = dict()
    for i, ec_number in enumerate(ec_numbers):
        print('\tMapping EC Number {} of {}'.format(i + 1, len(ec_numbers)))
        try:
            ec_response = json.loads(requests.get(constants.UNIPROT_EC_ENDPOINT.format(taxcode, ec_number)).text).get(
                'results')
            ec_to_uniprot[ec_number] = [entry.get('primaryAccession') for entry in ec_response]
        except json.decoder.JSONDecodeError:
            print('Error occurred with EC Number: {}'.format(ec_number))
    return ec_to_uniprot


def compute_kegg_node_label(node: etree.Element,
                            graphics: etree.Element,
                            compounds_drugs_and_glycans_to_name_dict: dict[str],
                            ko_to_name_dict: dict[str]) -> str | None:
    if node.get('type') == 'compound':
        return "; ".join(
            [compounds_drugs_and_glycans_to_name_dict[name.split(':')[-1]] for name in node.get('name').split() if
             name.split(':')[-1] in compounds_drugs_and_glycans_to_name_dict])
    elif node.get('type') == 'ortholog':
        return "; ".join([ko_to_name_dict[name] for name in node.get('name').split() if name in ko_to_name_dict])
    elif node.get('type') == 'map':
        return graphics.get('name')
    else:
        return None


def convert_kegg_kgml_to_json(kgml_file: Path, compounds_drugs_and_glycans_to_name_dict: dict[str],
                              ko_to_name_dict: dict[str],
                              all_kegg_gene_ids: set) -> dict:
    """
    A KEGG KGML file consists in 'entry', 'reaction' and 'relation' tags (as well as a toplevel 'pathway' tag containing metadata).
    An 'entry' can be a simple node, e.g. gene, compound, or map (link to another pathway), or it can be a 'group'.
    A 'reaction' and a 'relation' both describe a link between two or more nodes, albeit in a slightly different notation.
    We first process the metadata, then all non-group 'entry' tags, then the group entries, and finally relations and reactions.

    Args:
        kgml_file: The raw XML definition of the KEGG pathway
        compounds_drugs_and_glycans_to_name_dict: A dictionary that maps KEGG compound identifiers to human-readable names
        ko_to_name_dict: A dictionary that maps KEGG ortholog identifiers to human-readable names
        all_kegg_gene_ids: A set of all appearing KEGG gene identifiers that is incrementally built while parsing all pathways.
                            This set is afterwards mapped to Uniprot identifiers

    Returns: A JSON object consisting of three entries:
        - pathway (Object, contains metadata on the pathway)
        - nodes (List)
        - links (List)

    """
    tree = etree.parse(kgml_file)
    root = tree.getroot()

    res = dict()

    res['pathway'] = {
        'name': root.get('name'),
        'org': root.get('org'),
        'number': root.get('number'),
        'title': root.get('title'),
        'image': root.get('image'),
        'link': root.get('link')
    }

    res['nodes'] = []
    for entry in root.findall('entry'):
        # Exclude those nodes which do not have x/y coordinates as well as group nodes
        # (they are dealt with in the next step)
        if (graphics := entry.find('graphics')).get('x') is not None and entry.get('type') != 'group':
            node_name = graphics.get('name')

            # Skip 'TITLE' nodes, we don't want them in the visualization
            # (it would just be the name of the pathway randomly floating around)
            if node_name.startswith('TITLE'):
                continue

            # Trim away ellipsis if present
            if node_name.endswith('...'):
                node_name = node_name[:-3]

            res['nodes'].append({
                'id': entry.get('id'),
                'defaultName': node_name.split(', ')[0] if entry.get('type') in ['gene', 'ortholog'] else None,
                'geneNames': node_name.split(', ') if entry.get('type') in ['gene', 'ortholog'] else [],
                'keggIds': entry.get('name').split(),
                'type': constants.NODE_TYPE_DICTIONARY.get(entry.get('type'), 'misc'),
                'link': entry.get('link'),
                'x': int(graphics.get('x')),
                # Experience has shown that the KEGG graphs look better when scaled by 1.5 vertically
                'y': int(graphics.get('y')) * 1.5,
                'label': compute_kegg_node_label(entry, graphics, compounds_drugs_and_glycans_to_name_dict,
                                                 ko_to_name_dict),
            })
            if entry.get('type') == 'gene':
                all_kegg_gene_ids |= {*entry.get('name').split()}

    # Now scan the group nodes
    group_entries = [entry for entry in root.findall('entry') if entry.get('type') == 'group']
    for entry in group_entries:
        components = [c.get('id') for c in entry.findall('component')]
        # Only add a group if it has members
        if components:
            res['nodes'].append({
                'id': entry.get('id'),
                # Should always be 'group' but just for consistency
                'type': constants.NODE_TYPE_DICTIONARY.get(entry.get('type'), 'misc'),
                # TODO: See if this comes in handy, if not can be deleted again
                'components': components
            })
        for c in entry.findall('component'):
            # ...add the group ids to every group member (to the 'res' dictionary)
            component_node = [node for node in res.get('nodes') if node.get('id') == c.get('id')][0]
            component_node['groupId'] = entry.get('id')

    # Some group nodes have shared members, which does not make sense (e.g. hsa05014, entries 449 and 451).
    # We merge them, meaning we only retain one of them, add all nodes from one into the other (if they are not already in there)
    # And change all edge source/targets from one to the other
    group_ids_to_replace = {}
    for group1 in res.get('nodes'):
        if group1.get('id') in group_ids_to_replace or group1.get('type') != 'group':
            continue
        for group2 in res.get('nodes'):
            if group1.get('id') == group2.get('id') or group2.get('id') in group_ids_to_replace or group2.get(
                    'type') != 'group':
                continue
            if len(set(group1.get('components')).intersection(group2.get('components'))) > 0:
                group1['components'] = list(set(group1['components']).union(group2.get('components')))
                group_ids_to_replace[group2.get('id')] = group1['id']

    # Drop all groups with an 'id to replace'
    res['nodes'] = [node for node in res['nodes'] if node['id'] not in group_ids_to_replace.keys()]

    # Update the groupIds of all nodes
    for node in res['nodes']:
        if 'groupId' in node:
            node['groupId'] = group_ids_to_replace.get(node['groupId'], node['groupId'])

    # Now process the links (relation/reaction tags).
    # We might have duplicates afterwards.
    # So we keep track of the already added
    # source-target combinations.
    res['links'] = []
    unique_links = set()
    for i, relation in enumerate(root.findall('relation')):
        # Check on the fly if it is a 'to-replace' group id, if not, just use it as is
        link_tuple = ((source_id := group_ids_to_replace.get(relation.get('entry1'), relation.get('entry1'))),
                      (target_id := group_ids_to_replace.get(relation.get('entry2'), relation.get('entry2'))))
        if link_tuple not in unique_links:
            unique_links.add(link_tuple)
            link_subtypes = [constants.LINK_TYPE_DICTIONARY.get(subtype.get('name'), subtype.get('name'))
                             for subtype in relation.findall('subtype')]
            link_labels = [constants.LINK_LABEL_DICTIONARY.get(subtype) for subtype in link_subtypes]
            res['links'].append({
                'id': 'relation-{}'.format(i),
                'sourceId': source_id,
                'targetId': target_id,
                'types': link_subtypes,
                'label': ", ".join([label for label in link_labels if label])
            })

    for i, reaction in enumerate(root.findall('reaction')):
        for j, substrate in enumerate(reaction.findall('substrate')):
            for k, product in enumerate(reaction.findall('product')):
                link_tuple = ((source_id := group_ids_to_replace.get(substrate.get('id'), substrate.get('id'))),
                              (target_id := group_ids_to_replace.get(product.get('id'), product.get('id'))))
                if link_tuple not in unique_links:
                    unique_links.add(link_tuple)
                    res['links'].append({
                        'id': 'reaction-{}_{}_{}'.format(i, j, k),
                        'name': reaction.get('name'),
                        'sourceId': source_id,
                        'targetId': target_id,
                        'type': 'reaction',
                        'types': [reaction.get('type')],
                    })

    return res


def convert_wikipathways_gpml_to_json(gpml_file: Path, gene_ids_by_database: defaultdict[set],
                                      organism_identifier: str) -> dict:
    """
    A Wikipathways GPML file consists in a multitude of different tags.
    Relevant for us are these types (in the order in which they are processed):
     - 'DataNode': Can e.g. be Genes, Proteins, Metabolites,
        contains 'Graphics' with the node's coordinates and 'Xref' with external database identifier(s)
     - 'Group': Consists of identifiers which can be referenced by DataNodes or Interactions
     - 'Interaction': Contains source and target identifiers (can be DataNodes, Groups, or even other Interactions)
            Also contains an 'ArrowHead' property which we use to classify the type of interaction.
            If an interaction has as source or target another interaction, this is realized using an 'Anchor' tag inside
            the Interaction.
     - 'Label': Any kind of text label. Some (e.g. links to other pathways) can be part of interactions or groups.
            In that case, they are added as Nodes to the JSON, otherwise they are skipped.
     - 'State': Indiciates e.g. Ubiquitylation of a protein. In the JSON, states will be appended to the label of the
            DataNode which they belong to.

    Args:
        gpml_file: The raw XML definition of the Wikipathways pathway
        gene_ids_by_database: A dictionary that maps all supported gene identifier databases to sets of all identifiers
            that appear in Wikipathways. This dictionary of sets is incrementally built while parsing all pathways.
            Each set is afterwards mapped to Uniprot identifiers.
        organism_identifier: Identifier (e.g. 'hsa', 'eco') of the organism currently being processed

    Returns: A JSON object consisting of three entries:
        - pathway (Object, contains metadata on the pathway)
        - nodes (List)
        - links (List)

    """

    # Some elements might be missing a graph id, so as a fallback, we keep our own counter
    fallback_id = 0

    ### 0. Preliminary processing
    tree = etree.parse(gpml_file)
    root = tree.getroot()

    res = dict()

    res['pathway'] = {
        'name': gpml_file.name.split('_')[-2],
        'org': root.get('Organism'),
        'title': root.get('Name'),
    }

    # The set of databases that appears in Wikipathways is quite heterogenous. We need to map some of the database names
    # to aliases that are understood by the Uniprot API.
    database_name_map = {
        'BRENDA': 'EC',
        'Entrez Gene': 'GeneID',
        'Enzyme Nomenclature': 'EC',
        'KEGG Genes': 'KEGG',
        'NCBI Protein': 'RefSeq_Protein',
        'RefSeq': 'RefSeq_Protein',
        'Uniprot-TrEMBL': 'Uniprot',
    }

    def get_database_name(db: str, id: str) -> str:
        """
        Usually, this is just a dictionary lookup.
        There are two exceptions, HGNC and Ensembl,
        which need to be sub-categorized based on the format of the id
        """
        if db in ['HGNC', 'HGNC Accession number']:
            return 'HGNC' if id.startswith('HGNC:') or id.isnumeric() else 'GeneCards'
        elif id.startswith('ENSG'):
            return 'Ensembl'
        elif id.startswith('ENST'):
            return 'Ensembl_Transcript'
        else:
            return database_name_map.get(db, db)

    def format_gene_id(id_raw: str, db: str, organism_identifier: str) -> Union[str, None]:
        """The gene identifiers are also not unified, so we need to look out for some edge cases."""
        if not id_raw:
            return None
        elif db == 'HGNC' and ':' not in id_raw:
            return 'HGNC:{}'.format(id_raw)
        elif db == 'KEGG':
            return '{}:{}'.format(organism_identifier, id_raw)
        else:
            return id_raw.strip()

    ### 1. DataNode tags ###
    # Some nodes belong to a group. The groups have two identifiers:
    # Like any other nodes, a GraphId by which they are referenced in Links
    # But they also have a groupId, by which they are referenced in their group members
    # We will get rid of the groupIds by mapping them to their graphIds
    groupId2GraphId = {groupnode.get('GroupId'): groupnode.get('GraphId')
                       for groupnode in root.xpath("//*[starts-with(local-name(), 'Group')]") if
                       groupnode.get('GraphId') and groupnode.get('GroupId')}
    nodesDict = {}
    for datanode in root.xpath("//*[starts-with(local-name(), 'DataNode')]"):
        if len(graphicss := datanode.xpath(".//*[starts-with(local-name(), 'Graphics')]")) > 0:
            references = datanode.xpath(".//*[starts-with(local-name(), 'Xref')]")

            # Collect the reference ID of the node so they can later all be mapped together
            reference_database_raw = references[0].get('Database') if len(references) > 0 else None
            reference_id_raw = references[0].get('ID') if len(references) > 0 else None
            reference_database = get_database_name(reference_database_raw, reference_id_raw)
            reference_id = format_gene_id(reference_id_raw, reference_database, organism_identifier)

            if reference_database in constants.GENE_ID_DATABASES and reference_id:
                gene_ids_by_database[reference_database].add(reference_id)

            if datanode.get('GraphId'):
                node_id = datanode.get('GraphId')
            else:
                node_id = str(fallback_id)
                fallback_id += 1

            nodesDict[node_id] = {
                'id': node_id,
                'geneNames': [datanode.get('TextLabel')] if constants.NODE_TYPE_DICTIONARY.get(
                    datanode.get('Type')) == constants.GENE_NODETYPE else None,
                'type': constants.NODE_TYPE_DICTIONARY.get(datanode.get('Type'), 'misc'),
                'x': float(graphicss[0].get('CenterX')),
                'y': float(graphicss[0].get('CenterY')),
                'label': datanode.get('TextLabel') if constants.NODE_TYPE_DICTIONARY.get(
                    datanode.get('Type')) != constants.GENE_NODETYPE else None,
                'groupId': groupId2GraphId.get(datanode.get('GroupRef')),
                'referenceDatabase': reference_database,
                'referenceId': reference_id,
            }

    ### 2. Interaction tags ###
    # We need to pass over the links another time and possibly delete some of them, so we start with a temporary list
    links_temp = []
    # Anchors, which are used in the GPML structure to realize Interactions between Interactions, have their own IDs.
    # We do not include these anchors in the JSON, instead we directly map Interactions to one another.
    # Therefore we need to create a mapping of anchor IDs to Interaction IDs
    anchor_graphId_to_interaction_graphId = dict()
    # Now iterate over the interactions
    for interaction in root.xpath("//*[starts-with(local-name(), 'Interaction')]"):
        if len(graphicss := interaction.xpath(".//*[starts-with(local-name(), 'Graphics')]")) > 0:
            graphics = graphicss[0]
            points = [point for point in graphics.xpath(".//*[starts-with(local-name(), 'Point')]")
                      if point.get('GraphRef') is not None]

            # Some edges have only one endpoint, we cannot interpret them and skip them for now
            if len(points) != 2:
                continue

            source, target = (points[0], points[1]) if points[1].get('ArrowHead') is not None else (
                points[1], points[0])
            links_temp.append({
                'id': interaction.get('GraphId'),
                'sourceId': source.get('GraphRef'),
                'targetId': target.get('GraphRef'),
                'types': [constants.LINK_TYPE_DICTIONARY.get(target.get('ArrowHead'))]
            })

            if not links_temp[-1]['id']:
                links_temp[-1]['id'] = fallback_id
                fallback_id += 1

            # If the interaction has an anchor, add its ID to the map of anchor to interaction
            for anchor in graphics.xpath(".//*[starts-with(local-name(), 'Anchor')]"):
                anchor_graphId_to_interaction_graphId[anchor.get('GraphId')] = interaction.get('GraphId')

    # Now we do a second pass over the interactions to replace the anchor ids by actual source/target ids
    # using the map we built in the first pass:
    res['links'] = []
    unique_links = set()
    for link in links_temp:
        # If an end of an edge is another edge, replace its ID (which points to the anchor) by the edge ID itself
        link['sourceId'] = anchor_graphId_to_interaction_graphId.get(link['sourceId'], link['sourceId'])
        link['targetId'] = anchor_graphId_to_interaction_graphId.get(link['targetId'], link['targetId'])
        # Check that the link is not a duplicate (could happen e.g. if both source and target were a group node before)
        link_tuple = (link["sourceId"], link["targetId"])
        if link_tuple not in unique_links:
            unique_links.add(link_tuple)
            res["links"].append(link)

    ### 3. Label tags ###
    for label in root.xpath("//*[starts-with(local-name(), 'Label')]"):
        addLabel = False
        for source, target in unique_links:
            if label.get('GraphId') == source or label.get('GraphId') == target:
                addLabel = True
                break
        # Patch: Skip labels that only have 'P' as text label, they would be confusing for the user
        if not addLabel and (groupOfLabel := label.get('GroupRef')) is not None and label.get('TextLabel') != 'P':
            for groupId in groupId2GraphId:
                if groupOfLabel == groupId:
                    addLabel = True
                    break
        if addLabel:
            if len(graphicss := label.xpath(".//*[starts-with(local-name(), 'Graphics')]")) > 0:
                if label.get('GraphId'):
                    node_id = label.get('GraphId')
                else:
                    node_id = str(fallback_id)
                    fallback_id += 1
                nodesDict[node_id] = {
                    'id': node_id,
                    'geneNames': None,
                    'type': constants.NODE_TYPE_DICTIONARY.get(label.get('Type'), 'misc'),
                    'x': float(graphicss[0].get('CenterX')),
                    'y': float(graphicss[0].get('CenterY')),
                    'label': label.get('TextLabel'),
                    'groupId': groupId2GraphId.get(label.get('GroupRef'))
                }

    ### 4. State tags ###
    # We add the state tags as a list to their associated Data Nodes
    # Right now, the biowcpathwaygraph is not using this information, but in the future it might.
    for state in root.xpath("//*[starts-with(local-name(), 'State')]"):
        if 'states' not in nodesDict[state.get('GraphRef')]:
            nodesDict[state.get('GraphRef')]['states'] = []
        nodesDict[state.get('GraphRef')]['states'].append(state.get('TextLabel'))

    ### 5. Group tags ###
    # Add each group as a node
    for groupnode in root.xpath("//*[starts-with(local-name(), 'Group')]"):
        if groupnode.get('GraphId') and groupnode.get('GroupId'):
            node_id = groupnode.get('GraphId')
            # In case a node has been assigned a fallback Id, we want it in here.
            components = [node.get('id') for node in nodesDict.values() if
                          node.get('groupId') == groupId2GraphId.get(groupnode.get('GroupId'))]
            # Only add a group if it has members
            if components:
                nodesDict[node_id] = {
                    'id': node_id,
                    'type': 'group',
                    'components': components
                }
    res['nodes'] = list(nodesDict.values())

    return res


def get_kegg_idlist_to_names(idtype: str) -> dict[str, str]:
    print('Retrieving list of KEGG {} names...'.format(idtype))
    api_response = requests.get(constants.KEGG_GENERIC_IDENTIFIER_LIST_ENDPOINT.format(idtype))
    identifier_to_name_df = pd.read_csv(BytesIO(api_response.content), sep='\t', names=['identifier', 'name'])
    return dict(zip(identifier_to_name_df['identifier'], identifier_to_name_df['name']))


def process_kegg_pathways(raw_dir: Path, json_dir: Path, organism_identifier: str) -> None:
    compounds_drugs_and_glycans_to_name_dict = get_kegg_idlist_to_names('cpd')
    compounds_drugs_and_glycans_to_name_dict.update(get_kegg_idlist_to_names('dr'))
    compounds_drugs_and_glycans_to_name_dict.update(get_kegg_idlist_to_names('gl'))
    ko_to_name_dict = get_kegg_idlist_to_names('ko')

    all_kegg_gene_ids = set()
    all_jsondata = dict()
    pathway_ids = get_kegg_pathway_list(organism_identifier)
    new_pathway_ids = []
    for pathway_id in tqdm(pathway_ids, 'Converting KEGG Pathways to JSON'):
        raw_kgml_outputfile = Path(raw_dir, '{}.kgml'.format(pathway_id))
        if not raw_kgml_outputfile.exists():
            new_pathway_ids.append(pathway_id)
            continue
        all_jsondata[pathway_id] = convert_kegg_kgml_to_json(raw_kgml_outputfile,
                                                             compounds_drugs_and_glycans_to_name_dict,
                                                             ko_to_name_dict,
                                                             all_kegg_gene_ids)

    print('Mapping KEGG Gene Identifiers to Uniprot...')
    kegg_to_uniprot = map_gene_identifiers('KEGG', 'UniProtKB', all_kegg_gene_ids)

    all_uniprot_ids = set()
    for pathway_name, pathway_data in tqdm(all_jsondata.items(), 'Adding Uniprot IDs'):
        for node in pathway_data['nodes']:
            if node['type'] == 'gene_protein':
                # Pop the keggIds because we don't need them anymore after this
                uniprot_ids = {uniprot_id for kegg_id in node.pop('keggIds') for uniprot_id in
                               kegg_to_uniprot.get(kegg_id, [])}
                node['uniprotAccs'] = list(uniprot_ids)
                all_uniprot_ids |= uniprot_ids

    print('Mapping Uniprot identifiers to Gene Names...')
    uniprot_to_gene_names = map_gene_identifiers('UniProtKB_AC-ID', 'Gene_Name', all_uniprot_ids)

    json_dir.mkdir(exist_ok=True, parents=True)
    for pathway_name, pathway_data in tqdm(all_jsondata.items(), 'Adding Gene Names and writing out JSON Files'):
        for node in pathway_data['nodes']:
            if node['type'] == 'gene_protein':
                node['geneNames'] = list(
                    set(node['geneNames']) | {gene_name for uniprot in node['uniprotAccs'] for gene_name in
                                              uniprot_to_gene_names[uniprot]})
        json_outputfile = Path(json_dir, '{}.json'.format(pathway_name))
        with open(json_outputfile, 'w') as o:
            json.dump(pathway_data, o, indent=2)

    if new_pathway_ids:
        print('New KEGG pathways found on the server. Consider running the script again with --download.\n'
              'These are the new pathway IDs: {}'.format(new_pathway_ids))


class OrganismInfo(TypedDict):
    identifier: str
    taxcode: int


def process_wikipathways_pathways(raw_dir: Path, json_dir: Path, organism: OrganismInfo) -> None:
    gene_ids_by_database = defaultdict(set)
    all_jsondata = dict()
    for gpml_file in tqdm(list(Path.glob(raw_dir, '*.gpml')), 'Converting Wikipathways Pathways to JSON'):
        pathway_name = gpml_file.stem
        all_jsondata[pathway_name] = convert_wikipathways_gpml_to_json(gpml_file, gene_ids_by_database,
                                                                       organism['identifier'])

    mapped_gene_ids_by_database = dict()
    for source_db, gene_ids in gene_ids_by_database.items():
        print('Wikipathways: Mapping {} identifiers to Uniprot'.format(source_db))
        if source_db == 'Uniprot':
            continue
        elif source_db == 'EC':
            mapped_gene_ids_by_database[source_db] = map_ec_to_uniprot(organism['taxcode'], gene_ids)
        # For some reason it doesn't work for Refseq with the full thing (the request never terminates but
        # does not give an error either)
        # So we split it into several parts and update the results dictionary; since the input is a set we do not
        # run a risk of overwriting stuff
        elif source_db == 'RefSeq_Protein':
            mapped_gene_ids_by_database[source_db] = defaultdict(list)
            parts = 10
            ids_mapped = 0
            for i in range(parts):
                print('RefSeq part {}/{}'.format(i + 1, parts))
                mapped_gene_ids_by_database[source_db].update(
                    map_gene_identifiers(source_db, 'UniProtKB', set(list(gene_ids)[i::parts])))
                ids_mapped += len(set(list(gene_ids)[i::parts]))
                print('Mapped {} of {} identifiers.'.format(ids_mapped, len(gene_ids)))
        else:
            mapped_gene_ids_by_database[source_db] = map_gene_identifiers(source_db, 'UniProtKB', gene_ids)

    json_dir.mkdir(exist_ok=True, parents=True)
    for pathway_name, pathway_data in tqdm(all_jsondata.items(), 'Adding Uniprot IDs and writing out JSON Files'):
        for node in pathway_data['nodes']:
            if node.get('referenceDatabase') in constants.GENE_ID_DATABASES:
                if node.get('referenceDatabase') == 'Uniprot':
                    node['uniprotAccs'] = [node['referenceId']] or []
                else:
                    node['uniprotAccs'] = mapped_gene_ids_by_database[node['referenceDatabase']].get(
                        node['referenceId'])

        wp_identifier = pathway_name.split('_')[-2]
        json_outputfile = Path(json_dir, '{}.json'.format(wp_identifier))
        with open(json_outputfile, "w") as o:
            json.dump(pathway_data, o, indent=2)

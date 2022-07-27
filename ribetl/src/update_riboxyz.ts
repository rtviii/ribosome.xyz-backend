// https://search.rcsb.org/#introduction
import axios from "axios";
import { gzip, ungzip } from 'node-gzip'
import fs from 'fs'
import yargs from 'yargs'
import path from "path";
import shell from "shelljs";
import { processPDBRecord } from "./requestGqlProfile";
import { RibosomeStructure } from "./RibosomeTypes";
import { BlobOptions } from "buffer";



const missing_structures = async () => {
    console.log("Getting missing structures");
    var rcsb_search_api = "https://search.rcsb.org/rcsbsearch/v2/query"
    const params = {
        "query":
        {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "group",
                    "logical_operator": "and",
                    "nodes": [
                        { "type": "group", "logical_operator": "and", "nodes": [{ "type": "terminal", "service": "text", "parameters": { "operator": "contains_phrase", "negation": false, "value": "RIBOSOME", "attribute": "struct_keywords.pdbx_keywords" } }] },
                        { "type": "group", "logical_operator": "and", "nodes": [{ "type": "terminal", "service": "text", "parameters": { "operator": "greater", "negation": false, "value": 25, "attribute": "rcsb_entry_info.polymer_entity_count_protein" } }] },
                        { "type": "group", "logical_operator": "and", "nodes": [{ "type": "terminal", "service": "text", "parameters": { "operator": "less", "negation": false, "value": 4, "attribute": "rcsb_entry_info.resolution_combined" } }] }
                    ],
                    "label": "text"
                }
            ],
            "label": "query-builder"
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": true,
            "results_verbosity": "compact"
        }
    };
    let query = rcsb_search_api + "?json=" + encodeURIComponent(JSON.stringify(params))

    // let cypherstring = "match (struct:RibosomeStructure) return struct.rcsb_id"
    // cypherstring = encodeURIComponent(cypherstring);
    // let ribxz_query = `http://localhost:8000/neo4j/cypher/?cypher=${cypherstring}`

    let dbquery = new Promise<string[]>((resolve, reject) => {
        shell.config.silent = true

        console.log("Executing :\n\t\t",`echo \"match (struct:RibosomeStructure) return struct.rcsb_id\" | cypher-shell -a \"${process.env.NEO4J_URI}\" --format plain -u ${process.env.NEO4J_USER} -p ${process.env.NEO4J_PASSWORD} --database ${process.env.NEO4J_CURENTDB}`)
        shell.exec(`echo \"match (struct:RibosomeStructure) return struct.rcsb_id\" | cypher-shell -a \"${process.env.NEO4J_URI}\" --format plain -u ${process.env.NEO4J_USER} -p ${process.env.NEO4J_PASSWORD} --database ${process.env.NEO4J_CURENTDB}`,
            function (err, stdout, stderr) {
                if (err != 0) {
                    console.log("Got shell error code: ", err)
                    reject(err)
                }

                const dbstructs = (stdout as string).replace(/"/g, '').split("\n").filter(r => r.length === 4)
                resolve(dbstructs)
            })
    })



    return Promise.all([dbquery, axios.get(query)]).then(r => {
        var ribxz_structs: string[] = r[0]
        var rcsb_structs : string[] = r[1].data.result_set
        

        var missing_from_ribxz = rcsb_structs.filter(struct => {
            if (!ribxz_structs.includes(struct)) {
                return true
            } else return false
        })
        console.log(`riboxyz contains ${ribxz_structs.length} structures. Up-to-date RCSB API contains ${rcsb_structs.length} structures.`)
        console.log("Structs absent from ribosome.xyz: ", missing_from_ribxz.length)
        return missing_from_ribxz
    })
    .catch(e => { console.log(`Got error: ${e}`); return [] })
}


const download_unpack_place = async (struct_id: string) => {
    const BASE_URL = "http://files.rcsb.org/download/"
    const FORMAT = ".cif.gz"

    const structid = struct_id.toUpperCase()
    let url = BASE_URL + structid + FORMAT
    let compressed: Buffer = await axios.get(url, { responseType: 'arraybuffer' }).then(r => { return r.data })
        .catch(e => { console.log(`Structure ${structid} failed: `, e); return []; })
    let decompressed = await ungzip(compressed);

    let destination_chains = path.join(
        process.env["STATIC_ROOT"] as string,
        `${structid}`,
        `CHAINS`)

    let structfile = path.join(
        process.env["STATIC_ROOT"] as string,
        `${structid}`,
        `${structid}.cif`)
    if (!fs.existsSync(destination_chains)) {
        fs.mkdirSync(destination_chains)
        console.log(`Created directory ${destination_chains}.`);
    }
    fs.writeFileSync(structfile, decompressed)
}

const save_struct_profile = (r: RibosomeStructure) => {
    var rcsb_id = r.rcsb_id;
    var target_filename = path.join(
        process.env.STATIC_ROOT as string,
        rcsb_id.toUpperCase(),
        rcsb_id.toUpperCase() + ".json"
    );

    if (!fs.existsSync(path.dirname(target_filename))) {
        shell.mkdir("-p", path.dirname(target_filename));
    }
    fs.writeFileSync(target_filename, JSON.stringify(r, null, 4));
    console.log("Has written to ", path.dirname(target_filename));
};




// Options describing how to process a given structure
interface IngressOptions{
    acquirePDBRecord   : boolean;
    downloadCifFile    : boolean;
    splitRenameChains  : boolean;
    extractBindingSites: boolean;
    commitToNeo4j      : boolean;
}

const process_new_structure = async (struct_id: string, envfilepath: string) => {
    struct_id = struct_id.toUpperCase()
    let ribosome = await processPDBRecord(struct_id)
    await save_struct_profile(ribosome)
    console.log(`Saved structure profile ${struct_id}.json`);
    await download_unpack_place(struct_id)
    console.log(`Saved structure cif file ${struct_id}.cif`);
    shell.exec(`${process.env.PYTHONPATH} /home/rxz/dev/riboxyzbackend/ribetl/src/split_rename.py -s ${struct_id} -env ${envfilepath}`)// renaming chains
    shell.exec(`${process.env.PYTHONPATH} /home/rxz/dev/riboxyzbackend/ribetl/src/bsite_mixed.py -s ${struct_id} -env ${envfilepath} --save`)// binding sites
    

    shell.exec(`export RIBOXYZ_DB_NAME="${process.env.NEO4J_DB_NAME}"; /home/rxz/dev/riboxyzbackend/ribetl/src/structure.sh ${path.join(process.env.STATIC_ROOT as string, struct_id, `${struct_id}.json`)}`)// binding sites
    shell.exec(`export RIBOXYZ_DB_NAME="${process.env.NEO4J_DB_NAME}"; /home/rxz/dev/riboxyzbackend/ribetl/src/proteins.sh ${path.join(process.env.STATIC_ROOT as string, struct_id, `${struct_id}.json`)}`)// binding sites
    shell.exec(`export RIBOXYZ_DB_NAME="${process.env.NEO4J_DB_NAME}"; /home/rxz/dev/riboxyzbackend/ribetl/src/rna.sh     ${path.join(process.env.STATIC_ROOT as string, struct_id, `${struct_id}.json`)}`)// binding sites
    shell.exec(`export RIBOXYZ_DB_NAME="${process.env.NEO4J_DB_NAME}"; /home/rxz/dev/riboxyzbackend/ribetl/src/ligands.sh ${path.join(process.env.STATIC_ROOT as string, struct_id, `${struct_id}.json`)}`)// binding sites
}

const main = async () => {
    // https://github.com/yargs/yargs/blob/main/docs/typescript.md
    const args = yargs(process.argv.slice(2)).options({
        envfile   : { type: 'string', demandOption: true },
        structure : { type: "string", demandOption: false, alias: "s", },
        pythonpath: { type: "string", demandOption: false, alias: "pypath", },

        neo4juri     : { type: "string", demandOption: false, alias: "neo4juri", },
        password: { type: "string", demandOption: false, alias: "p", },
        user    : { type: "string", demandOption: false, alias: "u", },
        dbname  : { type: "string", demandOption: false, alias: "d", },

        dryrun: { type: "boolean", demandOption: false, alias: "dryrun", },
    })
        .boolean('all')
        .parseSync();

    const DEFAULT_PYTHON_PATH = "/home/rxz/dev/riboxyzbackend/venv/bin/python3"
    process.env.PYTHONPATH = args.pythonpath || DEFAULT_PYTHON_PATH

    require('dotenv').config({ path: args['envfile'] });

    if (args.neo4juri ) { process.env.NEO4J_URI      = args.neo4juri }
    if (args.dbname   ) { process.env.NEO4J_CURENTDB = args.dbname   }
    if (args.user     ) { process.env.NEO4J_USER     = args.user     }
    if (args.password ) { process.env.NEO4J_PASSWORD = args.password }



    if (args.dryrun) {
        console.log("Dry run. No changes will be made to the database.")
        let missing = await missing_structures()
        // console.log("Missing structures: ", missing);
        // process.exit(1)
        missing.forEach(m => console.log(m))


    }

    if (args.all) {
        let missing = await missing_structures()
        console.log(`Attempting to process ${missing.length} structures`)
        missing.forEach(async (struct_id) => {
            await process_new_structure(struct_id, args.envfile)
        })
        process.exit(1)
    }


    if (args.structure) {
        if (args.structure.length < 4) {
            console.log("Enter a valid RCSB PDB ID to build from RCSB's gql api.");
            process.exit(2);
        }
        console.log(`Processing ${args.structure}`);
        try {
            await process_new_structure(args.structure, args.envfile);
        } catch (e: any) {
            console.log("Failed: \n\n");
            console.log(e);
        }
    }
}

main()
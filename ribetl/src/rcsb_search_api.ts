// https://search.rcsb.org/#introduction
import axios from "axios";
import { gzip, ungzip } from 'node-gzip'
import fs from 'fs'
var http = require('http');
import { createGzip, unzip } from 'zlib'






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
    let cypherstring = "match (struct:RibosomeStructure) return struct.rcsb_id"
    cypherstring = encodeURIComponent(cypherstring);

    let ribxz_query = `http://localhost:8000/neo4j/cypher/?cypher=${cypherstring}`

    return Promise.all([axios.get(ribxz_query), axios.get(query)]).then(r => {
        var ribxz_structs: string[] = r[0].data
        var rcsb_structs: string[] = r[1].data.result_set

        var missing_from_ribxz = rcsb_structs.filter(struct => {
            if (!ribxz_structs.includes(struct)) {
                return true
            } else return false
        })
        console.log(`riboxyz contains ${ribxz_structs.length} structures. Up-to-date RCSB API contains ${rcsb_structs.length} structures.`)
        console.log("Structs absent from ribosome.xyz: ", missing_from_ribxz.length)
        return missing_from_ribxz
    }).catch(e => { console.log(`Rejected : ${e}`); return [] })
}

const place_unpack_downloads = async (rcsb_structs_to_download: string[]) => {


    const BASE_URL = "http://files.rcsb.org/download/"
    const FORMAT = ".cif.gz"

    for (var structid of rcsb_structs_to_download.slice(0, 5)) {
        let url = BASE_URL + structid.toUpperCase() + FORMAT
        let compressed:Buffer = await axios.get(url, {responseType:'arraybuffer'}).then(r => {return r.data})
            .catch(e => {console.log(`Structure ${structid} failed: `, e); return [];})

        let decompressed = await  ungzip(compressed);
        fs.writeFileSync(`./${structid}.cif`, decompressed)

        // let _ = unzip(compressed, (err: any, data) => { console.log("decompressed data", data); fs.writeFileSync(`./${structid}.cif`, data.toString()) });


    }
    // let gzips = await Promise.all(pending_downloads);
    // console.log(`Got gzips : ${gzips}`)
    // let uncompressed = await Promise.all(gzips.map(async gz => await ungzip(gz)) as Promise<Buffer>[])

    // uncompressed.forEach( (v)=>{

    // }
    // )



    // return fs.writeFile(`./${structid}.cif`, await ungzip(r.data), function (err) {
    //     if (err) return console.log(err);
    //     console.log(`written to file ${structid}`);
    // });
}



const update_ribosome_xyz = async () => {
    await place_unpack_downloads(await missing_structures())
}


update_ribosome_xyz().then()
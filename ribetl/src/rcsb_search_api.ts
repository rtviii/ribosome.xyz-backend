// https://search.rcsb.org/#introduction
import axios from "axios";

var rcsb_search_api = "https://search.rcsb.org/rcsbsearch/v2/query"
const params = { "query": 
{   "type": "group",
    "logical_operator": "and", 
    "nodes": [
        {   "type": "group",
            "logical_operator": "and",
            "nodes": [
                      { "type": "group", "logical_operator": "and", "nodes": [{ "type": "terminal", "service": "text", "parameters": { "operator": "contains_phrase", "negation": false, "value": "RIBOSOME", "attribute": "struct_keywords.pdbx_keywords" } }]},
                      { "type": "group", "logical_operator": "and", "nodes": [{ "type": "terminal", "service": "text", "parameters": { "operator": "greater", "negation": false, "value": 25, "attribute": "rcsb_entry_info.polymer_entity_count_protein" } }] },
                      { "type": "group", "logical_operator": "and", "nodes": [{ "type": "terminal", "service": "text", "parameters": { "operator": "less", "negation": false, "value": 4, "attribute": "rcsb_entry_info.resolution_combined" } }] }
                    ], 
            "label": "text" }
            ], 
"label"        :   "query-builder" },
"return_type"  :   "entry"        ,
"request_options": {
    "return_all_hits": true,
    "results_verbosity": "compact"
}
};
let query = rcsb_search_api + "?json=" + encodeURIComponent(JSON.stringify(params))
axios.get(query).then(r => console.log(r.data.result_set)).catch(e => console.log(e))



let our_structs = 




// axios.get(search_api, {
//     params,
//     paramsSerializer: (p: any) => qs.stringify(p, { encode: true })
// }).then(r => console.log(r.data)).catch(e => { console.log(`Could fetch it. Err: ${e}`) })

// axios.get("https://search.rcsb.org/rcsbsearch/v2/query?json=%7B%22query%22%3A%7B%22type%22%3A%22group%22%2C%22logical_operator%22%3A%22and%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22full_text%22%2C%22parameters%22%3A%7B%22value%22%3A%22%5C%22thymidine%20kinase%5C%22%22%7D%7D%2C%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22operator%22%3A%22exact_match%22%2C%22value%22%3A%22Viruses%22%2C%22attribute%22%3A%22rcsb_entity_source_organism.taxonomy_lineage.name%22%7D%7D%2C%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22operator%22%3A%22exact_match%22%2C%22value%22%3A%22X-RAY%20DIFFRACTION%22%2C%22attribute%22%3A%22exptl.method%22%7D%7D%2C%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22operator%22%3A%22less_or_equal%22%2C%22value%22%3A2.5%2C%22attribute%22%3A%22rcsb_entry_info.resolution_combined%22%7D%7D%2C%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22operator%22%3A%22greater%22%2C%22attribute%22%3A%22rcsb_entry_info.nonpolymer_entity_count%22%2C%22value%22%3A0%7D%7D%5D%7D%2C%22return_type%22%3A%22entry%22%7D").then(r=>console.log(r.data)).catch(e=>{console.log(`Could fetch it. Err: ${e}`)})


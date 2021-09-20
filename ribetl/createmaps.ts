import {large_subunit_map} from './resources/large-subunit-map'
import * as fs from 'fs'
import {small_subunit_map} from './resources/small-subunit-map'



fs.writeFileSync("SSUMap.json", JSON.stringify(small_subunit_map, null, 4))
fs.writeFileSync("LSUMap.json", JSON.stringify(large_subunit_map, null, 4))
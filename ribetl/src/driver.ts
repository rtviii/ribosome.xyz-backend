import yargs from "yargs";
import shell from "shelljs";
import { processPDBRecord } from "./requestGqlProfile";
import { RibosomeStructure } from "./RibosomeTypes";
import path from "path";
import fs from 'fs'

export const writeupdateStruct = (r: RibosomeStructure) => {
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


(async () => {
})();



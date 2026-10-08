// Detail rendering + locus graphic. Depends on globals from index.html.
const UCSC={hs1:"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hs1&position=",
  hg38:"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=",
            t2t:"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hs1&position="};
/* Reverse the build's shard packing (build_dashboard.pack_repeats, layout 1) into
   the row shapes the renderers use. Idempotent; needs LOOKUP loaded. */
function unpackLocus(d){
  if(!d||d._unpacked) return d;
  const LK=(typeof LOOKUP!=="undefined"&&LOOKUP)||{};
  if(d.rep&&!d.repeats){
    const names=LK.repnames||[], rows=[];
    for(const asm of Object.keys(d.rep)){
      const c=d.rep[asm]; let p=0;
      for(let i=0;i<c.s.length;i++){
        p+=c.s[i]; const t=names[c.n[i]]||["?","?","?"];
        rows.push({assembly:asm,start:p,end:p+c.l[i],strand:c.d[i],rep_name:t[0],rep_class:t[1],
                   rep_family:t[2],pct_div:c.v[i]/10,n_loci:c.k[i]});
      }
    }
    d.repeats=rows;
  }
  if(!d.group_info) d.group_info=(LK.groups||{})[d.group]||{};
  if(!d.superfamily_info) d.superfamily_info=(LK.superfamilies||{})[d.group_info.superfamily]||{};
  d._unpacked=true;
  return d;
}
function coordOf(d,asm){return (d.coord||[]).find(c=>c.assembly===asm);}

// ---- v0.2 transcriptional units -------------------------------------------
// The TU layer is a SEPARATE assignment from the locus layer: units are built
// from element walks, not from loci, so one locus can span several units and a
// unit can carry several loci. Absent for loci in no released v0.2 unit.
function tuPanel(d){
  const t=d.tu||{};
  if(!t.units||!t.units.length){
    const why = t.in_tu===false ? "this locus is not represented in the v0.2 unit set"
                                 : "no v0.2 unit assignment in this bundle";
    return '<div class="panel"><h2>Transcriptional unit (v0.2)</h2>'+
           '<div class="note">'+esc(why)+'</div></div>';
  }
  const frag=[];
  if(t.n_tu_overlap_hg38>1) frag.push("hg38: spans "+t.n_tu_overlap_hg38+" units"+
    (t.dom_frac_hg38!=null?", dominant covers "+(100*t.dom_frac_hg38).toFixed(0)+"% of the locus":""));
  if(t.n_tu_overlap_t2t>1)  frag.push("t2t: spans "+t.n_tu_overlap_t2t+" units"+
    (t.dom_frac_t2t!=null?", dominant covers "+(100*t.dom_frac_t2t).toFixed(0)+"% of the locus":""));
  const warn = frag.length
    ? '<div class="note">'+esc(frag.join(" \u00b7 "))+
      ' \u2014 the locus is not coextensive with any single unit; the dominant unit is listed first</div>'
    : '';
  const rows=t.units.map(function(u){
    return '<tr><td><span class="cid" style="font-size:13px">'+esc(u.tu_id)+'</span></td>'+
      '<td>'+esc(u.group_call||"")+'<span class="note"> ('+esc(u.call_level||"")+')</span></td>'+
      '<td>'+esc(u.verdict||"")+'</td>'+
      '<td>'+esc(u.provenance||"")+'</td>'+
      '<td>'+esc(u.assemblies||"")+'</td>'+
      '<td>'+(u.internal_bp_hg38!=null?(+u.internal_bp_hg38).toLocaleString():
              (u.internal_bp_t2t!=null?(+u.internal_bp_t2t).toLocaleString()+' <span class="note">(t2t)</span>':""))+'</td>'+
      '<td>'+esc(u.evidence||"")+'</td>'+
      '<td>'+(u.is_chimeric?'<b>chimeric</b> ':'')+esc(u.member_groups||"")+'</td></tr>';
  }).join("");
  return '<div class="panel"><h2>Transcriptional unit (v0.2)</h2>'+
    '<div class="note">A separate assignment from the locus layer: units are built from '+
    'element walks, so a locus may span several units and a unit may carry several loci. '+
    'Cite the unit id only alongside the locus versioned_id.</div>'+warn+
    '<table class="tbl" style="margin-top:8px"><thead><tr><th>unit</th><th>group call</th>'+
    '<th>verdict</th><th>provenance</th><th>assemblies</th><th>internal bp</th>'+
    '<th>evidence</th><th>member groups</th></tr></thead><tbody>'+rows+'</tbody></table></div>';
}

function render(d,h){
  if(!d){$("view").innerHTML='<div class="empty">locus not found in bundle</div>';return;}
  d=unpackLocus(d);
  if(d.uid&&typeof history!=="undefined"&&history.replaceState&&typeof location!=="undefined"&&location.hash!=="#"+d.uid)
    try{ history.replaceState(null,"","#"+d.uid); }catch(_){}
  if(typeof LASTUID!=="undefined"&&d.uid){ LASTUID=d.uid; if(typeof setMode==="function") setMode("locus"); }
  const hg=coordOf(d,"hg38"),t2=coordOf(d,"t2t"),g=d.group_info||{},st=d.structure||{},xg=d.crossgenome||{};
  // graphic assembly: hg38 when available, else t2t. The t2t fallback is what gives the
  // ~2,945 loci absent from hg38 a graphic at all; their lanes are repeats-only.
  const gco=hg||t2, gasm=hg?"hg38":(t2?"t2t":null);
  const links=[hg?'<a class="ucsc" target="_blank" href="'+UCSC.hg38+hg.chrom+":"+(hg.start+1-1000)+"-"+(hg.end+1000)+'">UCSC hg38 ±1 kb</a>':"",
               t2?'<a class="ucsc" target="_blank" href="'+UCSC.t2t+t2.chrom+":"+(t2.start+1-1000)+"-"+(t2.end+1000)+'">UCSC T2T (hs1) ±1 kb</a>':""].join("");
  const dfb=(d.dfam_best||[])[0]||{};
  $("view").innerHTML=
   '<div class="panel"><div class="idline"><span class="cid">'+esc(d.combined_id)+'</span>'+
     '<span class="uid">'+esc(d.uid)+'</span><span class="uid">'+esc(d.versioned_id)+'</span></div>'+
     '<div class="note">Cite the versioned_id. locus_uid is the immutable primary key; combined_id is positional and may re-letter.</div>'+
     '<div style="margin-top:9px">'+links+fastaButtons(d,hg,t2)+'</div>'+
     '<div class="note" id="fastanote" style="min-height:0"></div></div>'+
   '<div class="panel"><h2>Locus map ±1 kb ('+(gasm||"—")+')'+
     (gco&&gasm!=="hg38"?' <span class="note" style="font-weight:400">— hg38 coordinate absent;'+
       ' gene models, gEVE ORFs and HERVarium domains are hg38-only and are omitted</span>':"")+
     (gco?'<span class="dlbtns"><button class="ucsc dl" data-dl="svg" title="download this map as SVG '+
       '(drawn in your browser)">SVG</button><button class="ucsc dl" data-dl="png" title="download this map as PNG, '+
       '2\u00d7 resolution (drawn in your browser)">PNG</button></span>':"")+
     '</h2>'+
     '<div id="gfx">'+
     (gco?'<div class="note">rendering…</div>'
        :'<div class="note">no hg38 or t2t coordinate — graphic unavailable</div>')+'</div>'+
     '</div>'+
   '<div class="two">'+
     '<div class="panel"><h2>Locus</h2><dl class="kv">'+
       kv("group",d.group)+kv("band",d.band)+kv("origin",d.origin)+
       kv("structure",st.structure)+kv("category",st.category)+
       kv("LTR names",st.ltr_names)+kv("internal names",st.int_names)+
       kv("tandem/nested",st.is_tandem_or_nested?"yes":"no")+
       kv("segments",d.segments.length+(d.segments.length?"":" (none stored)"))+
       spliceKv(d)+
     '</dl></div>'+
     '<div class="panel"><h2>Group — '+esc(d.group)+'</h2><dl class="kv">'+
       kv("superfamily",g.superfamily)+kv("HERV class",g.herv_class)+
       kv("loci in group",g.n_loci)+kv("int model",g.intModel)+
       kv("RepBase class",g.repbase_class)+kv("HERVd family",g.hervd_family)+
       kv("dominant LTR",g.dominant_ltr)+
       kv("with flanking LTR",g.frac_with_flanking_ltr==null?null:(100*g.frac_with_flanking_ltr).toFixed(1)+"%")+
       kv("extension verdict",g.extension_verdict)+
       // Lineage is many-to-many with group, so it gets two rows, never one:
       // the viewed locus's own Navigator lineage(s) -- exact -- and the
       // group's dominant lineage carried WITH its share, so a reader can see
       // when "dominant" means 85% and when it means 38%.
       lineageKv(d,g)+
     '</dl></div></div>'+
   '<div class="panel"><h2>Aliases — '+new Set(d.aliases.map(a=>a.alias_type)).size+' types</h2>'+aliasTable(d.aliases)+'</div>'+
   tuPanel(d)+
   cclePanel(d)+
   tissuePanel(d)+rnaAtlasPanel(d)+
   // transcription evidence sits next to the TU panel: both describe what is
   // transcribed here, and the split question reads them together
   (typeof txPanel==="function"?txPanel(d):"")+
   (typeof jxPanel==="function"?jxPanel(d):"")+
   '<div class="panel"><h2>Coordinates</h2>'+coordTable(d,xg)+'</div>'+
   hml2Panel(d.hml2_detail)+
   missillacPanel(d.missillac)+
   (dfb.consensus_name?'<div class="panel"><h2>Dfam best alignment</h2><dl class="kv">'+
       kv("consensus",dfb.consensus_name)+kv("accession",dfb.dfam_accession)+
       kv("% identity",dfb.pct_identity==null?null:Number(dfb.pct_identity).toFixed(1))+
       kv("consensus coverage",dfb.cons_cov==null?null:(100*dfb.cons_cov).toFixed(1)+"%")+
       kv("SW score",dfb.sw_score)+kv("quality",dfb.aln_quality)+'</dl></div>':"")+
   (d.geve.length?'<div class="panel"><h2>gEVE ORFs — '+d.geve.length+'</h2>'+tbl(d.geve,
      ["geve_orf_id","orf_class","hg38_chrom","hg38_orf_start","hg38_orf_end","orf_strand"])+'</div>':"")+
   (d.domains.length?'<div class="panel"><h2>HERVarium domains — '+d.domains.length+'</h2>'+tbl(d.domains,
      ["gene","domain_desc","element","status","domain_score","hg38_start","hg38_end","strand"])+'</div>':"")+
   hervariumPanel(d)+hvorfPanel(d)+geneExprPanel(d)+
   (d.genes.length?'<div class="panel"><h2>Overlapping genes — '+d.genes.length+' rows</h2>'+tbl(d.genes,
      ["genome","source","ref_gene_name","ref_gene_id","gene_type","overlap_type","overlap_bp","exon_overlap_bp","n_transcripts"])+'</div>':"")+
   (d.segments.length?'<div class="panel"><h2>Segments — '+d.segments.length+'</h2>'+tbl(d.segments,
      ["seg_index","segment_class","repName","repFamily","repClass","chrom","start","end","strand","span","rmsk_sw_score"])+'</div>':"")+
   absentPanel(d,dfb);
  if(gco) drawLocus(d,gco,gasm);
  document.querySelectorAll("button.fa").forEach(b=>b.onclick=()=>fetchFasta(d,b.dataset.fa,$("fastanote")));
  document.querySelectorAll("button.dl").forEach(b=>b.onclick=()=>downloadMap(d,gasm,b.dataset.dl));
  document.querySelectorAll("button.tsv").forEach(b=>b.onclick=()=>{
    TS_ALL=b.dataset.ts==="1"; const c=$("tschart"); if(c) c.innerHTML=tsChart(d);
    document.querySelectorAll("button.tsv").forEach(x=>x.classList.toggle("on",x.dataset.ts===b.dataset.ts)); });
}
function kv(k,v){return "<dt>"+esc(k)+"</dt><dd>"+fmt(v)+"</dd>";}
/* CCLE segment shading: log2(fraction + 0.01), rescaled to 0-1, so the low end
   (1 vs 30 of 1,019 lines) is as distinguishable as the high end. 0 lines = outline. */
function ccO(n,NL){
  if(!n) return 0;
  const g=(Math.log2(n/NL+0.01)-Math.log2(0.01))/(Math.log2(1.01)-Math.log2(0.01));
  return 0.12+0.88*Math.max(0,Math.min(1,g));
}
/* CCLE expression panel: two independent readouts of the same 1,019 public CCLE
   RNA-seq runs. Telescope = EM read reassignment (TPM, hg38, Telescope loci only);
   31-mer = presence of locus-unique k-mers in each run (METHODS S17-S20). */
/* Per-tissue bar chart: share of each tissue's CCLE lines meeting a criterion.
   Shared by the Telescope (>= 1 TPM) and 31-mer (body detected) readouts so the
   two use identical axes and encoding. k[i] = lines meeting it in tissue i. */
function tisChart(TIS,k,extra,crit){
  const ML=30, BW=19, CH=60, W=ML+TIS.length*BW+80, H=CH+104, Y=f=>CH-f*CH+4;
  const short={haematopoietic_and_lymphoid_tissue:"blood / lymphoid",central_nervous_system:"CNS",
    upper_aerodigestive_tract:"upper aerodigestive"};
  const axis=[[0,"0"],[0.5,"50%"],[1,"100%"]].map(p=>'<line x1="'+ML+'" y1="'+Y(p[0])+'" x2="'+(ML+TIS.length*BW)+'" y2="'+Y(p[0])+
    '" stroke="var(--axis)" stroke-width="'+(p[0]===0?0.9:0.6)+'"'+(p[0]===0.5?' stroke-dasharray="4,3"':p[0]===1?' stroke-dasharray="1,3"':"")+'/>'+
    '<text x="'+(ML-4)+'" y="'+(Y(p[0])+3)+'" font-size="8" fill="var(--mut)" text-anchor="end">'+p[1]+"</text>").join("");
  const bars=TIS.map((T,i)=>{
    const n=T[1], c=k[i]||0, f=n?c/n:0, h=Math.max(f>0?1.5:0,f*CH), X=ML+i*BW, few=n<5;
    const nm=short[T[0]]||T[0].replace(/_/g," ");
    const rows=[["cell lines",n+(few?" (too few for a reliable share)":"")],[crit,c+" ("+(100*f).toFixed(0)+"%)"]]
      .concat((extra?extra(i):[]).filter(r=>r[1]!=null));
    return '<g'+tipA(card(T[0].replace(/_/g," "),rows,"var(--ccle)"))+'><rect x="'+X+'" y="0" width="'+BW+'" height="'+(CH+10)+'" fill="transparent"/>'+
      '<rect x="'+(X+2)+'" y="'+(Y(f)-(f>0&&f*CH<1.5?1.5-f*CH:0))+'" width="'+(BW-4)+'" height="'+h+'" fill="var(--ccle)" rx="1"'+
      (few?' fill-opacity="0.35"':"")+'/>'+
      '<text transform="translate('+(X+BW/2-2)+','+(CH+12)+') rotate(55)" font-size="8" fill="var(--mut)">'+esc(nm)+" ("+n+")</text></g>";
  }).join("");
  return '<svg class="tischart" viewBox="0 0 '+W+' '+H+'" width="100%" style="max-width:'+W+'px">'+axis+bars+"</svg>";
}
/* Primary-tissue panel (tissue_layer.py, METHODS S22). Per group: share of samples with
   the body detected in recount3 coverage (pale bar) and share with >= 3 own SENSE split
   reads in the body from Snaptron (dark bar; junction strand comes from the splice motif,
   so it is the element's own transcription even though these libraries are unstranded).
   Default view: priority comparisons + the locus's strongest groups; toggle: all groups. */
const TS_PRIORITY=[["blood",["GTEx Whole Blood","TCGA LAML tumor","CCLE haematopoietic & lymphoid"]],
  ["lymphoid",["GTEx Spleen","GTEx Cells - EBV-transformed lymphocytes","TCGA DLBC tumor"]],
  ["testis",["GTEx Testis","TCGA TGCT tumor"]]];
const TS_KIND={"GTEx normal":["GTEx normal","var(--cat1)"],"TCGA tumour":["TCGA tumour","var(--cat2)"],
  "TCGA adjacent normal":["TCGA adjacent normal","var(--cat6)"],"TCGA metastatic":["TCGA metastatic","var(--cat4)"],
  "cell line (CCLE)":["CCLE cell lines","var(--ccle)"]};
const TS_TIER={1:["A","own sense splicing up","var(--cat2)"],2:["B","excess over local background up","var(--cat5)"],
  3:["C","host gene / host splicing explains it","var(--cat6)"],4:["D","unresolved","var(--neutral)"]};
let TS_ALL=false;
function tsRows(d,idx,meta,head){
  const ts=d.ts, G=meta.groups, RH=23, LW=210, BW=260, HX=LW+BW+10, HW=44, C1=HX+HW+156, C2=C1+68, C3=C2+76, W=C3+300;
  // value columns (D46): median element RPKM, body exon : other body coverage, highest samples
  const co38=(d.coord||[]).find(c=>c.assembly==="hg38"), hasCov=!!(d.cov&&d.cv&&d.pf&&co38);
  const TTL=((LOOKUP||{}).top_samples||{}).labels||[], TT={}; (d.tt||[]).forEach(e=>{TT[e[0]]=e[1];});
  const fv=v=>v>=100?Math.round(v):v>=10?v.toFixed(0):v.toFixed(1);
  const C=meta.cx||{}, cx=ts.cx, lxOff=!!(cx&&C.primary&&C.primary[cx[0]]==="intronic"&&C.host_orient[cx[1]]==="antisense");
  const rk=v=>v==null?"\u2014":(v/100)<0.1?(v>0?"<0.1":"0"):(v/100).toFixed((v/100)<10?1:0);
  const pc=v=>v==null||v<0?"\u2014":v===0?"0%":v<10?"<1%":Math.round(v/10)+"%";
  const flt=["no unique flank segment \u2014 dominance cannot reject read-through here","upstream flank only","downstream flank only","both flanks"][ts.fl||0];
  let y=14, out='<text x="'+(HX+HW+6)+'" y="9" font-size="8.5" fill="var(--mut)">det \u00b7 dom \u00b7 local \u00b7 spliced</text>'+
    '<text x="'+C1+'" y="9" font-size="8.5" fill="var(--mut)">median RPKM</text>'+
    '<text x="'+C2+'" y="9" font-size="8.5" fill="var(--mut)">exon : rest</text>'+
    '<text x="'+C3+'" y="9" font-size="8.5" fill="var(--mut)">highest samples (mean coverage / base)</text>';
  for(const blk of idx){
    if(blk[0]){ out+='<text x="0" y="'+(y+11)+'" font-size="10" font-weight="600" fill="var(--ink)">'+esc(blk[0])+"</text>"; y+=16; }
    for(const i of blk[1]){
      const g=G[i], dv=ts.d[i], bv=ts.b?ts.b[i]:null, sv=ts.s[i], few=g[3]<10, k=TS_KIND[g[2]]||[g[2],"var(--neutral)"];
      const lv=(ts.lx&&!lxOff)?ts.lx[i]:null, hv=ts.he?ts.he[i]:null;
      const tip=card(g[0],[["kind",k[0]],["samples (coverage)",g[3]+(few?" \u2014 too few for a reliable share":"")],
        ["detected (cpb > 1)",pc(dv)],["body-dominant (> 2\u00d7 flanks)",pc(bv)+" \u00b7 "+flt],["samples (junctions)",g[4]||"no junction data"],
        ["\u2265 "+meta.splice_min_reads+" own sense split reads"+(meta.splice_excl?" (unique junctions)":""),pc(sv)],
        ...(ts.sa?[["\u2003including non-unique junctions",pc(ts.sa[i])]]:[]),
        ...(ts.sb?[["\u2003also excluding possibly shared",pc(ts.sb[i])]]:[]),
        ["local excess (\u2265 2\u00d7 flanking background)",lxOff?"not shown: intronic, host antisense (ENCODE: excess here is host-strand)":ts.lx?pc(lv):"no background window"],
        ["median element / host (RPKM)",rk(ts.em?ts.em[i]:null)+" / "+(ts.hm?rk(ts.hm[i]):"no host gene")],
        ["host expressed (\u2265 1 RPKM)",ts.he?pc(hv):"no host gene"]],k[1]);
      const bar=(v,yy,h,op)=>v>0?'<rect x="'+LW+'" y="'+(y+yy)+'" width="'+Math.max(1.5,BW*v/1000)+'" height="'+h+'" fill="'+k[1]+'" fill-opacity="'+(few?op*0.5:op)+'"/>':"";
      out+='<g'+tipA(tip)+'><rect x="0" y="'+y+'" width="'+W+'" height="'+RH+'" fill="transparent"/>'+
        '<text x="'+(LW-6)+'" y="'+(y+12)+'" font-size="9.5" text-anchor="end" fill="'+(few?"var(--mut)":"var(--ink)")+'">'+
          esc(g[1])+' <tspan fill="var(--mut)">('+g[3]+")</tspan></text>"+
        '<rect x="'+LW+'" y="'+(y+2)+'" width="'+BW+'" height="'+(RH-4)+'" fill="var(--maptrack)"/>'+
        bar(dv,2,4,0.3)+bar(bv,7,4,0.62)+(lv!=null?(lv>0?'<rect x="'+LW+'" y="'+(y+12)+'" width="'+Math.max(1.5,BW*lv/1000)+'" height="4" fill="none" stroke="'+k[1]+'" stroke-width="1"/>':""):"")+bar(sv,17,4,1)+
        (hv!=null?'<rect x="'+HX+'" y="'+(y+5)+'" width="'+HW+'" height="'+(RH-10)+'" fill="var(--maptrack)"/>'+(hv>0?'<rect x="'+HX+'" y="'+(y+5)+'" width="'+Math.max(1.5,HW*hv/1000)+'" height="'+(RH-10)+'" fill="var(--neutral)"/>':""):"")+
        '<text x="'+(HX+HW+6)+'" y="'+(y+14)+'" font-size="9" fill="var(--mut)">'+pc(dv)+" \u00b7 "+pc(bv)+" \u00b7 "+(lv!=null?pc(lv):"\u2014")+(sv>=0?" \u00b7 "+pc(sv):"")+"</text>"+
        '<text x="'+C1+'" y="'+(y+14)+'" font-size="9" fill="var(--ink)">'+rk(ts.em?ts.em[i]:null)+"</text>"+
        '<text x="'+C2+'" y="'+(y+14)+'" font-size="9" fill="var(--ink)">'+(()=>{ const er=hasCov?covExonRest(d,co38,i):null;
          return er?(er[1]>0?"\u00d7"+fv(er[0]/er[1]):(er[0]>0?"rest < .01":"\u2014")):"\u2014"; })()+"</text>"+
        '<text x="'+C3+'" y="'+(y+14)+'" font-size="8.5" fill="var(--ink)">'+(TT[i]?(()=>{ const a=TT[i], o=[];
          for(let q=0;q<a.length;q+=2) o.push(esc(TTL[a[q]])+' <tspan fill="var(--mut)">'+covF(a[q+1]/1000)+"</tspan>");
          return o.join(" \u00b7 "); })():'<tspan fill="var(--mut)">none \u2265 '+(((LOOKUP||{}).top_samples||{}).min||0.3)+"</tspan>")+"</text></g>";
      y+=RH;
    }
    y+=6;
  }
  const ax=[0,0.5,1].map(f=>'<line x1="'+(LW+BW*f)+'" y1="0" x2="'+(LW+BW*f)+'" y2="'+(y-6)+'" stroke="var(--axis)" stroke-width="0.5"'+
    (f===0.5?' stroke-dasharray="3,3"':"")+'/><text x="'+(LW+BW*f)+'" y="'+(y+6)+'" font-size="8.5" text-anchor="middle" fill="var(--mut)">'+
    (f*100)+"%</text>").join("");
  const hax=ts.he?'<text x="'+(HX+HW/2)+'" y="'+(y+6)+'" font-size="8.5" text-anchor="middle" fill="var(--mut)">host</text>':"";
  return '<svg class="tsbars" viewBox="0 0 '+W+" "+(y+10)+'" width="100%" style="max-width:'+W+'px">'+out+ax+hax+"</svg>";
}
function tsKey(){
  const sw=(op,lab)=>'<span style="display:inline-block;width:22px;height:6px;background:var(--cat1);opacity:'+op+';margin:0 4px 1px 10px;vertical-align:middle"></span>'+lab;
  const ol='<span style="display:inline-block;width:22px;height:5px;border:1px solid var(--cat1);margin:0 4px 1px 10px;vertical-align:middle"></span>local excess';
  const hs='<span style="display:inline-block;width:22px;height:8px;background:var(--neutral);margin:0 4px 1px 10px;vertical-align:middle"></span>host gene expressed';
  return '<div class="note" style="font-style:normal;margin:2px 0 6px">'+sw(0.3,"detected")+sw(0.62,"body-dominant")+ol+sw(1,"own sense splicing")+hs+
    ' <span style="color:var(--mut)">\u00b7 numbers: detected \u00b7 body-dominant \u00b7 local excess \u00b7 spliced</span>'+
    '<br><span style="color:var(--mut)"><b>median RPKM</b>: median element expression in the group (unique segments). '+
    '<b>exon : rest</b>: mean coverage of the locus body\u2019s 100-bp bins that overlap a GENCODE v50 exon over the other body bins (all reads; '+
    '\u2014 = no exon in the body). <b>highest samples</b>: the group\u2019s top three samples by mean coverage per unique body base per 10\u2079 aligned '+
    'bases, listed from 0.3.</span></div>';
}
function tsJunctions(d,meta){
  const ja=d.ts.ja; if(!ja||!ja.length||!meta.splice_excl) return "";
  const C=meta.splice_excl.codes, nm=u=>(typeof IDX!=="undefined"&&IDX&&IDX.uid2cid&&IDX.uid2cid[u])||u;   // index not loaded yet: show the uid
  const COL={1:"var(--cat2)",2:"var(--cat2)",3:"var(--cat5)",4:"var(--neutral)"};
  const rows=ja.filter(r=>r[0]>0).map(r=>{
    const part=r[4].length?r[4].map(u=>'<a href="#'+esc(u)+'">'+esc(nm(u))+"</a>").join(", "):"\u2014";
    return '<tr><td style="font-family:ui-monospace,Menlo,monospace;font-size:11.5px">'+esc(r[1])+"</td><td>"+r[2].toLocaleString()+
      " ("+(r[3]<10?"<1":Math.round(r[3]/10))+'%)</td><td><b style="color:'+COL[r[0]]+'">'+esc(C[r[0]])+"</b></td><td>"+part+"</td></tr>";
  }).join("");
  const more=ja.find(r=>r[0]===0), any=ja.some(r=>r[0]<=2&&r[0]>0);
  return '<details class="tsja"'+(any?" open":"")+'><summary><b>Non-unique own-strand junctions</b> \u2014 '+
    (any?"excluded from the splicing share":"kept, flagged")+"</summary>"+
    '<table class="t tsjat"><tr><th>junction (hg38)</th><th>split reads (share of locus)</th><th>flag</th><th>also at</th></tr>'+rows+"</table>"+
    (more?'<div class="note">+ '+more[2]+" more flagged junctions.</div>":"")+
    '<div class="note">Junctions whose 50 bases around the splice site match another junction within 2 mismatches. '+
    "<b>Shared between loci</b>: identical sequence, and the reads split between the positions in the same ratio in every tissue, as the aligner would split one read pool. "+
    "<b>Same reads counted at another locus</b>: identical per-sample read counts at both. "+
    "Both are left out of the splicing share; for reads counted twice inside this locus (tandem copies) one copy is kept. "+
    "<b>Possibly shared</b> (kept): identical sequence with a tissue-dependent split, or 1\u20132 mismatches with a stable split.</div></details>";
}
function tsHost(d,meta){
  const hg=d.ts.hg; if(!hg) return "";
  const r=hg[2], hi=r!=null&&r>=100*((meta.host||{}).rho_flag||0.5);
  const nm=esc(hg[1])+(hg[1]!==hg[0].split(".")[0]?' <span style="color:var(--mut)">('+esc(hg[0].split(".")[0])+")</span>":"");
  const rt=r==null?"not defined (no variation)":"\u03c1 = "+(r/100).toFixed(2);
  return '<div class="'+(hi?"tsflag":"note")+'" style="'+(hi?"":"font-style:normal;margin:2px 0 6px")+'">'+(hi?"<b>Element coverage follows its host gene.</b> ":"")+
    "Host gene (GENCODE v50, exons outside HERVs): "+nm+". Element vs host across samples, ranked within each group: "+rt+"."+
    (hi?" Detection here largely tracks host transcription; body dominance and own splicing are the element-specific readouts.":"")+"</div>";
}
function tsChart(d){
  const meta=(LOOKUP||{}).tissue, ts=d.ts, G=meta.groups, gi={}; G.forEach((g,i)=>gi[g[0]]=i);
  if(!TS_ALL){
    const used=new Set(), blocks=TS_PRIORITY.map(p=>[p[0],p[1].filter(n=>n in gi).map(n=>(used.add(gi[n]),gi[n]))]);
    const top=G.map((g,i)=>[ts.d[i],i]).filter(r=>!used.has(r[1])&&G[r[1]][3]>=10&&r[0]>0).sort((a,b)=>b[0]-a[0]).slice(0,6).map(r=>r[1]);
    if(top.length) blocks.push(["highest for this locus (groups of \u2265 10)",top]);
    return tsRows(d,blocks,meta);
  }
  const order=["GTEx normal","TCGA tumour","TCGA adjacent normal","TCGA metastatic","cell line (CCLE)"];
  const blocks=order.map(k=>[(TS_KIND[k]||[k])[0],G.map((g,i)=>[g,i]).filter(r=>r[0][2]===k).sort((a,b)=>a[0][1].localeCompare(b[0][1])).map(r=>r[1])]).filter(b=>b[1].length);
  return '<div class="tsall">'+tsRows(d,blocks.slice(0,1),meta)+tsRows(d,blocks.slice(1),meta)+"</div>";
}
function tsContext(d,meta){
  const ts=d.ts, C=meta.cx||{}, cx=ts.cx, en=ts.en, n=meta.encode_experiments;
  const [o,sf]=en||[0,null], fl=ts.fl||0;
  const genes=cx&&cx[6]?cx[6].split("|").map(g=>esc(g)).join(", "):"";
  const enc=en?" ENCODE stranded RNA-seq ("+n+" tissue samples): "+(sf!=null?sf+"% of the signal over the element is on its own strand":"no signal over the element")+".":"";
  const flank=fl===0?" This locus has no unique flank segment, so body dominance cannot reject read-through here.":"";
  const box=(title,body)=>'<div class="tsflag"><b>'+title+"</b> "+body+"</div>";
  if(!cx) return en?'<div class="note">'+enc+"</div>":"";
  const prim=C.primary[cx[0]], hor=C.host_orient[cx[1]], xcl=C.exonic_class[cx[2]], sup=C.support[cx[3]], gcl=C.gclass[cx[4]], own=cx[7]===1;
  const READS="Unstranded detection (pale bars) counts every read over the element\u2019s unique sequence, so it records host transcription here as well. "+
    "Body dominance and own sense splicing are the more element-specific readouts.";
  if(prim==="intronic"){
    const hs={same:"on the element\u2019s strand",antisense:"on the opposite strand",both:"on both strands"}[hor];
    let why=hor==="antisense"&&sf!=null&&sf<20?" ENCODE confirms the host strand: only "+sf+"% of the stranded signal here is the element\u2019s own.":
            hor==="same"?" Because the host is transcribed on the element\u2019s own strand, stranded data cannot separate the two here either.":
            enc;
    return box("Intronic locus \u2014 detection is not element-specific.","It lies in an intron of "+(genes||"a host gene")+" ("+esc(gcl)+"), transcribed "+hs+
      ", with no exon of any GENCODE v50 transcript over it. Detection at intronic loci follows each sample\u2019s intronic RNA content (host pre-mRNA, retained introns) "+
      "and varies across tissues for that reason. "+READS+why+flank);
  }
  if(prim==="exonic"){
    const sup_={canonical:"a canonical (MANE / Ensembl) transcript",basic:"a basic transcript",long_read:"long-read (TAGENE) models only",other:"other transcript models"}[sup]||"";
    const what={ltr_driven_chimeric:"a transcript starts in this element\u2019s LTR and continues into exons outside it (LTR-driven, chimeric)",
      ltr_driven_contained:"a transcript starts in this element\u2019s LTR and stays within it (LTR-driven, contained)",
      internal_initiated_chimeric:"a transcript starts inside this element (not in an LTR) and continues into exons outside it",
      internal_initiated_contained:"a transcript starts inside this element and stays within it",
      exonised:"the element lies in an exon of a transcript that starts elsewhere",
      antisense_only:"only antisense transcripts have exons over the element"}[xcl];
    const tail=" ("+esc(gcl)+(genes?": "+genes:"")+(sup_?"; "+sup_:"")+")."+(cx[5]?" It is part of a protein-coding gene\u2019s MANE Select exon.":"");
    if(own) return '<div class="note" style="margin:4px 0 8px">GENCODE v50 annotates this element as its own gene model'+tail.replace(/\.$/,"")+
      ". Detection here is the element\u2019s transcription."+enc+flank+"</div>";
    if(xcl==="exonised"||xcl==="antisense_only")
      return box("Exonic in a "+(xcl==="exonised"?"host":"antisense")+" transcript \u2014 detection may be host signal.","GENCODE v50: "+what+tail+" "+READS+enc+flank);
    return '<div class="note" style="margin:4px 0 8px">GENCODE v50: '+what+tail+enc+flank+"</div>";
  }
  const on={1:"flanking transcription runs antisense to the element",2:"flanking transcription is on the element\u2019s strand",3:"flanking transcription is on both strands",0:"no flanking transcription"}[o];
  return '<div class="note" style="margin:4px 0 8px">Intergenic: no GENCODE v50 transcript spans the element.'+(en?" In ENCODE stranded RNA-seq, "+on+
    (sf!=null?"; "+sf+"% of the signal over the element is on its own strand.":"."):"")+flank+"</div>";
}
function tsTvn(d,meta){
  const tv=d.ts.tv; if(!tv||!tv.length) return "";
  const rows=tv.slice().sort((a,b)=>b[1]-a[1]).map(t=>{
    const p=meta.projects[t[0]], call=t[2]===1?"higher in tumour":t[2]===-1?"lower in tumour":"\u2014", tr=TS_TIER[t[3]];
    const w=Math.min(60,Math.abs(t[1])/10*20), col=tr?tr[2]:"var(--neutral)";
    const bar='<svg width="130" height="10"><line x1="65" y1="0" x2="65" y2="10" stroke="var(--axis)" stroke-width="0.6"/>'+
      '<rect x="'+(t[1]<0?65-w:65)+'" y="2" width="'+Math.max(1,w)+'" height="6" fill="'+col+'"'+(t[2]===0?' fill-opacity="0.35"':"")+"/></svg>";
    return "<tr><td>"+esc(p[0])+' <span class="note" style="font-style:normal">('+p[1]+")</span></td><td>"+(t[1]>0?"+":"")+(t[1]/10).toFixed(1)+
      "</td><td>"+bar+"</td><td>"+call+"</td><td>"+(tr?'<b style="color:'+tr[2]+'">'+tr[0]+"</b> "+esc(tr[1]):"")+"</td></tr>";
  }).join("");
  return '<h3 style="margin-top:12px">Tumour vs matched adjacent normal (TCGA)</h3><table class="t tstvn"><tr><th>project (pairs)</th>'+
    "<th>log\u2082 change</th><th></th><th>call (q &lt; 0.01, \u2265 2-fold)</th><th>evidence for the change</th></tr>"+rows+"</table>"+
    '<div class="note">Mean paired change in body coverage, unadjusted. Evidence: <b>A</b> the element\u2019s own sense splicing also rises; '+
    "<b>B</b> its excess over flanking exon-free sequence rises (not used for intronic elements antisense to the host, where ENCODE "+
    "shows that excess is host-strand); <b>C</b> the host gene or junctions spanning the element rise by at least half as much; "+
    "<b>D</b> none of these. Tumours carry more intronic RNA than their normals in most projects (intron retention), so a gain without A or B "+
    "may be host transcription.</div>";
}
function tissuePanel(d){
  const meta=(LOOKUP||{}).tissue; if(!meta||!d.ts) return "";
  return '<div class="panel" id="tspanel"><h2>Primary tissues \u2014 TCGA, GTEx and CCLE (recount3, Snaptron)'+
    '<span class="dlbtns"><button class="ucsc tsv'+(TS_ALL?"":" on")+'" data-ts="0">priority groups</button>'+
    '<button class="ucsc tsv'+(TS_ALL?" on":"")+'" data-ts="1">all '+meta.groups.length+" groups</button></span></h2>"+
    tsContext(d,meta)+tsHost(d,meta)+tsKey()+'<div id="tschart">'+tsChart(d)+"</div>"+tsJunctions(d,meta)+
    '<div class="note">'+"Bars per group: share of samples with the body <b>detected</b> (recount3 coverage over the locus\u2019s unique segments, cpb > 1); share <b>body-dominant</b> (detected, and mean coverage over the body more than twice that over each unique flank, which rejects read-through from surrounding transcription); share with \u2265 "+meta.splice_min_reads+" split reads on the element\u2019s own strand with an end in its body "+
    "(Snaptron; summed over its junctions; junctions that cannot be assigned to one locus are left out \u2014 see the list below the chart). <b>Local excess</b> (outlined): detected, and body coverage at least twice that of exon-free sequence 2\u201310 kb either side; not shown for intronic elements antisense to their host, where ENCODE puts that excess on the host strand. <b>Host</b> (grey, right): share of samples expressing the host gene (GENCODE v50 exons outside HERVs, \u2265 1 RPKM). Numbers: detected \u00b7 body-dominant \u00b7 local excess \u00b7 spliced. Samples in brackets; grey labels have fewer than 10. "+
    meta.excluded_samples+" samples without a tissue label or with very low coverage are excluded. "+
    "Detection in intron-rich samples partly reflects host pre-mRNA; own splicing does not.</div>"+
    tsTvn(d,meta)+"</div>";
}
function cclePanel(d){
  const cc=d.cc||{}, meta=(LOOKUP||{}).ccle; if(!meta||(!cc.tel&&!cc.bf)) return "";
  const NL=meta.n_lines, NK=meta.n_lines_kmer||NL,
    pctN=(n,D)=>n==null?"\u2014":n.toLocaleString()+" ("+(100*n/D).toFixed(n&&100*n/D<1?1:0)+"%)",
    pct=n=>pctN(n,NL), pctK=n=>pctN(n,NK);
  let tel="", telChart='<div></div>', kmChart='<div></div>';
  if(cc.tel){
    const t=cc.tel, TIS=meta.tissues;
    const chart=tisChart(TIS,t.tis.map(v=>v[2]),i=>[["at \u2265 1 TPM",null],["median TPM",t.tis[i][0]],["max TPM",t.tis[i][1]]],"at \u2265 1 TPM");
    tel='<div><h3>Telescope (TPM, hg38)</h3><div class="note" style="margin:-2px 0 4px">RetroelementDB quantification (Russ &amp; Iordanskiy, <i>Mobile DNA</i> 2025)</div><dl class="kv">'+kv("feature",t.id)+kv("max TPM",t.mx)+kv("median TPM",t.md)+
      kv("lines \u2265 0.1 TPM",pct(t.n01))+kv("lines \u2265 1 TPM",pct(t.n1))+kv("lines \u2265 10 TPM",pct(t.n10))+"</dl>"+
      (t.top.length?'<table class="cctop">'+t.top.map(r=>"<tr><td>"+esc(r[0])+'</td><td class="note" style="font-style:normal">'+
        esc(r[1].replace(/_/g," "))+"</td><td>"+r[2]+" TPM</td></tr>").join("")+"</table>":"")+
      "</div>";
    telChart='<div><div class="note" style="margin-top:2px">Telescope: share of each tissue\u2019s cell lines at \u2265 1 TPM '+
      '(number of lines in brackets; pale = fewer than 5)</div>'+chart+"</div>";
  } else tel='<div><h3>Telescope (TPM, hg38)</h3><div class="note">not a Telescope locus \u2014 Telescope\u2019s '+
    "annotation covers 14,968 HERV loci; this one is outside it</div></div>";
  let km="";
  const bf=cc.bf;
  if(bf&&bf.n_body!=null){
    const f5=bf.oriented?"5\u2032":"left", f3=bf.oriented?"3\u2032":"right";
    km='<div><h3>31-mer (locus-unique k-mers, '+esc(cc.asm)+')</h3><dl class="kv">'+
      kv("body detected in",pctK(bf.n_body))+
      kv(f5+" flank detected in",bf.n_5==null?"no unique flank segment":pctK(bf.n_5))+
      kv(f3+" flank detected in",bf.n_3==null?"no unique flank segment":pctK(bf.n_3))+
      kv("body > 2\u00d7 "+f5+" flank",bf.dom_5==null?"\u2014":pctK(bf.dom_5))+
      kv("body > 2\u00d7 "+f3+" flank",bf.dom_3==null?"\u2014":pctK(bf.dom_3))+
      kv("body > 2\u00d7 both flanks",bf.dom_both==null?"\u2014":pctK(bf.dom_both))+
      (bf.hs!=null?kv("high-specificity set","yes \u2014 spec_ds "+bf.hs):"")+"</dl>"+
      '<div class="note">Detected = the region\u2019s weighted sum of unique 31-mer counts exceeds '+meta.det_wsum+
      " (counts below "+meta.min_abund+" ignored; weight \u221d segment length). Against Telescope this call is "+
      "99.8% specific (lines Telescope scores 0) and 98% sensitive (lines at \u2265 1 TPM). "+'Body vs flank compares length-normalised weighted sums (weight \u221d segment length), so a long body '+
      "is not favoured over a 1 kb flank by size."+(bf.oriented?"":" No locus strand: flanks are genomic left/right.")+
      " Segment positions are on the map (CCLE 31-mer lane).</div>"+
      "</div>";
    if(bf.tis) kmChart='<div><div class="note" style="margin-top:2px">31-mer: share of each tissue\u2019s cell lines with the body '+
      "detected (weighted sum > "+meta.det_wsum+", as above)</div>"+tisChart(meta.tissues_kmer||meta.tissues,bf.tis,null,"body detected")+"</div>";
  } else km='<div><h3>31-mer (locus-unique k-mers)</h3><div class="note">no locus-unique 31-mer in the body: '+
    "k-mer evidence cannot be attributed to this locus (see mappability)</div></div>";
  return '<div class="panel"><h2>CCLE expression \u2014 '+NL.toLocaleString()+' cancer cell lines</h2><div class="ccgrid">'+
    tel+km+(cc.tel||(bf&&bf.tis)?telChart+kmChart:"")+'</div><div class="note" style="margin-top:8px">Two independent readouts of the same public CCLE RNA-seq runs: '+
    "Telescope reassigns multi-mapping reads by EM; the 31-mer readout asks whether k-mers unique to this locus occur in "+
    "each run. Where they disagree, the number of unique k-mers the locus has (mappability lanes, CCLE 31-mer lane) "+
    "is usually why."+(NK<NL?" The 31-mer readout uses "+NK.toLocaleString()+" of the "+NL.toLocaleString()+" runs: "+
      (NL-NK)+" runs whose k-mer records are ~25\u201330\u00d7 sparser than the rest (normal aligned coverage) are left out of it.":"")+
    "</div></div>";
}
/* FASTA of the locus body, fetched on click from the UCSC Genome Browser REST API.
   One request per click (guarded to at most one per second), strand-aware via
   revComp=1, coordinates passed as stored (0-based start, end-exclusive), which is
   the API's own convention. Nothing is prefetched. */
const UCSC_API="https://api.genome.ucsc.edu/getData/sequence";
const UCSC_DB={hg38:"hg38",t2t:"hs1"};
function fastaButtons(d,hg,t2){
  const b=(c,asm,lab)=>c?'<button class="ucsc fa" data-fa="'+asm+'" title="Fetches this locus\u2019s sequence from the '+
    'UCSC Genome Browser REST API (api.genome.ucsc.edu) when clicked; one request per click">FASTA '+lab+'</button>':"";
  return b(hg,"hg38","hg38")+b(t2,"t2t","T2T");
}
let _faLast=0;
async function fetchFasta(d,asm,note){
  const co=coordOf(d,asm); if(!co) return;
  const now=Date.now();
  if(now-_faLast<1000){ note.textContent="one request per second \u2014 try again in a moment"; return; }
  _faLast=now;
  const minus=co.strand==="-", unk=!(co.strand==="+"||co.strand==="-");
  const url=UCSC_API+"?genome="+UCSC_DB[asm]+";chrom="+encodeURIComponent(co.chrom)+
    ";start="+co.start+";end="+co.end+(minus?";revComp=1":"");
  note.textContent="fetching "+asm+" sequence from the UCSC Genome Browser REST API\u2026";
  let j;
  try{ const r=await fetch(url);
       if(!r.ok) throw new Error("HTTP "+r.status+(r.status===429?" (UCSC rate limit \u2014 wait and retry)":""));
       j=await r.json(); if(!j.dna) throw new Error(j.error||"no sequence in response"); }
  catch(e){ note.innerHTML="UCSC request failed: "+esc(e.message||String(e))+
    ' \u00b7 <a target="_blank" href="'+esc(url)+'">open the API URL</a>'; return; }
  const seq=j.dna.toUpperCase(), n=seq.length;
  if(n!==co.end-co.start){ note.textContent="UCSC returned "+n+" bp, expected "+(co.end-co.start)+" \u2014 not saved"; return; }
  const head=">"+d.versioned_id+" "+d.uid+" "+asm+" "+co.chrom+":"+(co.start+1)+"-"+co.end+
    "("+(unk?"+, locus strand unknown":co.strand)+") source=UCSC:"+UCSC_DB[asm];
  const fa=head+"\n"+seq.replace(/(.{60})/g,"$1\n").replace(/\n$/,"")+"\n";
  saveBlob(new Blob([fa],{type:"text/plain"}),d.combined_id+"_"+asm+".fa");
  note.innerHTML="saved "+esc(d.combined_id+"_"+asm+".fa")+" \u2014 "+n.toLocaleString()+" bp"+
    (minus?", reverse-complemented to the locus strand":"")+(unk?", + strand (locus strand unknown)":"")+
    ', from the <a target="_blank" href="https://genome.ucsc.edu/goldenPath/help/api.html">UCSC REST API</a>. '+
    "For many loci, use the bundle BEDs with the UCSC 2bit files instead.";
}
/* Shared download helper (tests replace window.__saveBlob). */
function saveBlob(blob,name){
  if(typeof window!=="undefined"&&window.__saveBlob) return window.__saveBlob(blob,name);
  const a=document.createElement("a"); a.href=URL.createObjectURL(blob); a.download=name;
  document.body.appendChild(a); a.click(); setTimeout(()=>{URL.revokeObjectURL(a.href); a.remove();},1000);
}
/* Standalone copy of the locus map: CSS variables resolved to literal colours (a
   downloaded SVG has no stylesheet), hover data stripped, white background, a title
   line and the legend drawn inside the image. Built entirely in the browser. */
function exportSVG(d,asm){
  const svg=document.querySelector("#gfx svg"); if(!svg) return null;
  const cs=getComputedStyle(document.documentElement);
  const rv=s=>s.replace(/var\(--([a-z0-9-]+)\)/g,(m,k)=>cs.getPropertyValue("--"+k).trim()||"#000");
  const c=svg.cloneNode(true);
  c.querySelectorAll("[data-tip]").forEach(e=>e.removeAttribute("data-tip"));
  c.querySelectorAll("title").forEach(e=>e.remove());
  c.querySelectorAll(".hit").forEach(e=>e.remove());
  c.querySelectorAll("*").forEach(e=>{ for(const a of ["fill","stroke","style"]){
    const v=e.getAttribute(a); if(v&&v.includes("var(")) e.setAttribute(a,rv(v)); } });
  const vb=(svg.getAttribute("viewBox")||"0 0 1080 200").split(/\s+/).map(Number);
  const W=vb[2], top=26, used=(document.getElementById("gfx").dataset.used||"").split(",").filter(Boolean);
  const LEG={ltr:"ERV LTR",int:"ERV internal",orf:"gEVE ORF",dom:"HERVarium domain",u3:"U3",rr:"R",u5:"U5",
    hvo:"HERVOminer ORF \u2265 81 aa (solid = gEVE)",gene:"gene exon / intron",
    cage:"FANTOM5 CAGE peak",["rep-line"]:"LINE",["rep-sine"]:"SINE",["rep-dna"]:"DNA",["rep-other"]:"other TE",
    ["rep-low"]:"simple / low complexity",map1:"uniquely mappable",ccle:"CCLE 31-mer segment"};
  let lx=110, ly=0, leg="";
  for(const k of Object.keys(LEG)){ if(!used.includes(k)) continue;
    const t=LEG[k], w=18+t.length*6.2;
    if(lx+w>W-10){ lx=110; ly+=16; }
    leg+='<rect x="'+lx+'" y="'+(ly+2)+'" width="10" height="10" rx="2" fill="'+rv("var(--"+k+")")+'"/>'+
      '<text x="'+(lx+14)+'" y="'+(ly+11)+'" font-size="10" fill="#444">'+esc(t)+"</text>"; lx+=w+10; }
  if(used.includes("jx")){ ly+=16;
    leg+='<text x="110" y="'+(ly+11)+'" font-size="10" fill="#444">junction samples</text>'+
      [[10,"10"],[1000,"1k"],[100000,"100k"]].map((p,i)=>'<line x1="'+(215+i*52)+'" y1="'+(ly+7)+'" x2="'+(235+i*52)+
      '" y2="'+(ly+7)+'" stroke="'+rv("var(--jx-within)")+'" stroke-width="'+jxW(p[0]).toFixed(2)+'" opacity="'+
      jxO(p[0]).toFixed(2)+'" stroke-linecap="round"/><text x="'+(239+i*52)+'" y="'+(ly+11)+'" font-size="9" fill="#666">'+
      p[1]+"</text>").join("")+'<text x="380" y="'+(ly+11)+'" font-size="10" fill="#444">above = sense, below = '+
      "antisense, dashed = no mappable anchor</text>"; }
  const H=vb[3]+top+ly+28, co=coordOf(d,asm)||{};
  const title=d.combined_id+"  ("+d.uid+", "+d.versioned_id+")  \u00b7  "+asm+" "+co.chrom+":"+
    ((co.start||0)+1).toLocaleString()+"-"+(co.end||0).toLocaleString()+
    (co.strand==="+"||co.strand==="-"?"  ("+(co.strand==="-"?"\u2212":"+")+" strand)":"");
  const inner=new XMLSerializer().serializeToString(c).replace(/^<svg[^>]*>/,"").replace(/<\/svg>\s*$/,"");
  return '<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" width="'+W+
    '" height="'+H+'" viewBox="0 0 '+W+" "+H+'" font-family="Menlo,Consolas,monospace">'+
    '<rect width="100%" height="100%" fill="#ffffff"/>'+
    '<text x="10" y="17" font-size="12" font-weight="bold" fill="#1b1b1f">'+esc(title)+"</text>"+
    '<g transform="translate(0,'+(top-vb[1])+')">'+inner+"</g>"+
    '<g transform="translate(0,'+(vb[3]+top+6)+')">'+leg+"</g>"+
    '<text x="'+(W-10)+'" y="'+(H-6)+'" font-size="9" fill="#888" text-anchor="end">HERV catalog v2 \u00b7 dashboard v1.0</text></svg>';
}
function downloadMap(d,asm,fmt){
  const s=exportSVG(d,asm); if(!s) return;
  const base=d.combined_id+"_"+asm+"_map";
  if(fmt==="svg"){ saveBlob(new Blob([s],{type:"image/svg+xml"}),base+".svg"); return; }
  const img=new Image(), url=URL.createObjectURL(new Blob([s],{type:"image/svg+xml"}));
  img.onload=()=>{ const k=2, cv=document.createElement("canvas");
    cv.width=img.width*k; cv.height=img.height*k;
    const g=cv.getContext("2d"); g.scale(k,k); g.drawImage(img,0,0); URL.revokeObjectURL(url);
    cv.toBlob(b=>saveBlob(b,base+".png"),"image/png"); };
  img.onerror=()=>{ URL.revokeObjectURL(url); alert("PNG export failed in this browser; use SVG."); };
  img.src=url;
}
/* One line naming the evidence this locus does not have, instead of a panel per
   absent layer. Each item says why it can be absent where that is structural. */
const LAYER_NAME={crossgenome:"hg38\u2013T2T alignment",dfam_aln:"Dfam alignment",dfam_hs1:"Dfam (hs1)",dfam_rm:"Dfam RepeatMasker",
  ervmap:"ERVmap",gene:"gene overlap",geve:"gEVE",hervarium:"HERVarium",hervd:"HERVd",repbase:"RepBase",retrotector:"RetroTector",
  tss:"TSS / FANTOM5",hervominer:"HERVOminer"};
function absentPanel(d,dfb){
  // "not screened" (the layer never looked at this locus) is kept apart from "none found"
  const ns=new Set(d.ns||[]), why=k=>ns.has(k)?"not screened":"none found";
  const miss=[];
  if(!dfb.consensus_name) miss.push(["Dfam alignment",ns.has("dfam_aln")?"not screened":"not screened against Dfam consensus"]);
  if(!d.geve.length) miss.push(["gEVE ORFs",why("geve")]);
  if(!d.domains.length) miss.push(["HERVarium domains",why("hervarium")]);
  if(!d.genes.length) miss.push(["overlapping genes",why("gene")]);
  if(!d.segments.length) miss.push(["stored segments","only Telescope-origin loci carry them"]);
  const nb=d.eb===2?'<div class="note" style="font-style:normal"><b>HERVarium-only locus</b> (entry batch 2): a HERVarium v5 element that '+
    "overlapped no earlier catalog locus (D41). Not screened by: "+[...ns].map(k=>LAYER_NAME[k]||k).join(", ")+
    "; nor by TU membership or the CCLE / TCGA / GTEx / ENCODE / Blueprint expression layers.</div>":"";
  if(!miss.length&&!nb) return "";
  return '<div class="panel"><h2>Not present for this locus</h2>'+nb+'<div class="absent">'+
    miss.map(m=>"<b>"+esc(m[0])+"</b>"+(m[1]?" ("+esc(m[1])+")":"")).join(" \u00b7 ")+"</div></div>";
}

// Snaptron splice-evidence summary. A strict ladder: the FIRST tier that
// applies is reported, so "donor or acceptor internal (sense)" implies no
// single sense junction had both ends inside, and an antisense tier implies no
// sense junction qualified at all.
//
// Two things this line deliberately does NOT do:
//  - it does not collapse "no evidence" into "not assessable". srav3h is
//    hg38-only and does not cover every alt contig, so ~4,000 loci cannot be
//    assessed at all; reporting those as negative would be a false statement.
//  - it does not derive the tier from the drawn arcs. Arcs are capped per
//    locus for bundle size, and the cap changes the tier on 405 loci, so the
//    line reports the full filtered evidence and says so when they can differ.
const SPLICE_TIER = {
  both_sense:  "donor and acceptor internal (sense)",
  one_sense:   "donor or acceptor internal (sense)",
  both_anti:   "donor and acceptor internal (antisense)",
  one_anti:    "donor or acceptor internal (antisense)",
  intronic:    "intronic \u2014 spanning junction"
};
const SPLICE_DETAIL = {
  donor:    "donor internal",
  acceptor: "acceptor internal",
  both:     "both, on separate splice junctions"
};
function spliceKv(d){
  const s=d.splice||{}, t=s.t;
  if(!t) return kv("splicing (Snaptron)",null);
  if(t==="not_assessable")
    return "<dt>splicing (Snaptron)</dt><dd><span class=\"note\">not assessable"+
           (s.d?" \u2014 "+esc(s.d):"")+"</span></dd>";
  if(t==="none")
    return "<dt>splicing (Snaptron)</dt><dd>no qualifying junction"+
           "<span class=\"note\"> (assessable; canonical, \u226510 samples)</span></dd>";
  let txt=SPLICE_TIER[t]||t;
  // The "specify which end" detail the ladder asks for, on the single-end tiers.
  if((t==="one_sense"||t==="one_anti") && s.d) txt+=" \u2014 "+(SPLICE_DETAIL[s.d]||s.d);
  // For the spanning tier the strand is what needs specifying, not the end.
  if(t==="intronic" && s.d) txt+=" ("+(s.d==="both"?"sense and antisense":esc(s.d))+")";
  return "<dt>splicing (Snaptron)</dt><dd>"+esc(txt)+
    (s.sc!=null?'<span class="note"> \u00b7 max '+(+s.sc).toLocaleString()+" samples</span>":"")+
    '<span class="note" title="Ladder: both-sense &gt; one-sense &gt; both-antisense &gt; '+
    'one-antisense &gt; intronic. Donor = the junction\u2019s 5\u2032 end. Computed from all '+
    'canonical junctions with \u226510 samples, not only the arcs drawn below \u2014 the '+
    'arc set is capped per locus.">\u00a0\u24d8</span></dd>';
}

function lineageKv(d,g){
  // This locus's own lineage(s), from the Navigator records overlapping it.
  const rows=(d.missillac||[]);
  const mine=[...new Set(rows.map(r=>r.lineage_id).filter(Boolean))];
  let out = mine.length
    ? "<dt>lineage (this locus)</dt><dd>"+mine.map(esc).join(", ")+
      (mine.length>1?'<span class="note"> \u00b7 '+mine.length+" Navigator records disagree</span>":"")+"</dd>"
    : kv("lineage (this locus)",null);
  if(g.dom_lineage){
    const pct=g.dom_frac==null?null:(100*g.dom_frac).toFixed(0);
    const weak=g.dom_frac!=null&&g.dom_frac<0.5;
    out += "<dt>group lineage</dt><dd>"+esc(g.dom_lineage)+
      '<span class="note'+(weak?" warn":"")+'"> \u00b7 '+
      (pct!=null?pct+"% of group loci":"dominant")+
      (g.n_lineages?", of "+g.n_lineages+" lineages in group":"")+"</span></dd>";
  }
  return out;
}
function tbl(rows,cols){
  if(!rows||!rows.length)return '<div class="note">none</div>';
  const use=cols.filter(c=>rows.some(r=>r[c]!=null&&r[c]!==""));
  return "<table><tr>"+use.map(c=>"<th>"+esc(c)+"</th>").join("")+"</tr>"+
    rows.map(r=>"<tr>"+use.map(c=>'<td class="'+(typeof r[c]==="number"?"mono":"")+'">'+fmt(r[c])+"</td>").join("")+"</tr>").join("")+"</table>";
}
/* Subramanian 2011 HML-2 provirus detail. Only 87 of 39,733 loci carry this, so
   the panel is omitted entirely rather than rendered empty. Coordinates shown are
   the paper's own hg19 -- the hg38/T2T equivalents are already in the Coordinates
   panel, and repeating them here would imply the paper published them. */
/* ERV Navigator (Missillac) panel.

   Several rows on one locus is normal, not an error: our spans are internal-element
   extents, so one longer Navigator element can contain two of our loci and one of
   our loci can be hit by an element plus its flanking LTR records. The primary row
   is the best Jaccard; `margin` exposes how thin that call was, because a 0.001
   margin between two candidates is a coin-flip dressed as a decision.

   The lineage reconciliation is shown per row with its purity, and a name/coordinate
   DISAGREEMENT is called out explicitly. Navigator names its lineages after the LTR
   family; this catalog groups by internal element, so e.g. lineage ERVR-1.Theta.MER21
   sits on loci we group as MER4B and none of them carry a MER21* RepBase name. That
   is a real vocabulary difference, not a mapping error, and hiding it would let a
   user read "MER21" as our MER21. */
function missillacPanel(rows){
  if(!rows||!rows.length)return "";
  const pri=rows.find(r=>r.is_primary)||rows[0];
  const navlink=r=>{
    const u=r.url||("https://ervnavigator.fredhutch.org/locus/"+encodeURIComponent(r.missillac_id));
    return '<a href="'+esc(u)+'" target="_blank" rel="noopener noreferrer">'+esc(r.missillac_id)+" \u2197</a>";
  };
  const disagree=pri.name_agrees_with_coords===0&&pri.name_match_level!=="numeric_expansion";
  const head='<dl class="kv">'+
    kv("category",pri.category)+
    kv("lineage",pri.lineage_id)+
    kv("clade",pri.clade)+
    kv("our group (by overlap)",pri.majority_group==null?null:
        pri.majority_group+" ("+(100*Number(pri.purity)).toFixed(0)+"% of "+
        (pri.confidence||"?")+"-confidence lineage hits)")+
    kv("lineage name matches",pri.name_match_group==null?"no name match":
        pri.name_match_group+" ("+(pri.name_match_level||"").replace(/_/g," ")+")")+
    kv("records on this locus",rows.length)+
    '</dl>';
  const warn=disagree
    ? '<div class="note" style="border-left:3px solid #a33;padding-left:8px">'+
      'Vocabulary conflict: the lineage name resolves to <b>'+esc(pri.name_match_group)+
      '</b> but the loci this lineage overlaps are grouped <b>'+esc(pri.majority_group)+
      '</b>. ERV Navigator names lineages after the LTR family; this catalog groups by '+
      'internal element. Treat the two names as different vocabularies, not synonyms.</div>'
    : "";
  const body=rows.map(r=>{
    const b=[];
    if(r.is_primary)b.push('<span class="badge">primary</span>');
    if(r.is_primary&&r.primary_margin!=null&&Number(r.primary_margin)<0.05)
      b.push('<span class="badge ambig" title="best and runner-up Jaccard differ by '+
        Number(r.primary_margin).toFixed(3)+' — a thin call">thin margin</span>');
    if(r.n_loci_for_record>1)b.push('<span class="badge ambig">'+r.n_loci_for_record+" loci</span>");
    if(r.strand_agree===0)b.push('<span class="badge collapse">strand differs</span>');
    if(r.lift_status&&r.lift_status!=="both")b.push('<span class="badge collapse">'+esc(r.lift_status)+"</span>");
    return "<tr><td>"+navlink(r)+"</td><td class=\"mono\">"+esc(r.rbrt_id||"")+"</td>"+
      "<td>"+esc(r.category||"")+"</td><td>"+esc(r.lineage_id||"")+"</td>"+
      '<td class="mono">'+(r.jaccard==null?"":Number(r.jaccard).toFixed(3))+"</td>"+
      '<td class="mono">'+(r.ovl_bp==null?"":Number(r.ovl_bp).toLocaleString())+"</td>"+
      "<td>"+b.join(" ")+"</td></tr>";
  }).join("");
  return '<div class="panel"><h2>ERV Navigator \u2014 '+rows.length+
    (rows.length===1?" record":" records")+"</h2>"+head+warn+
    "<table><tr><th>Missillac ID</th><th>RBRT</th><th>category</th><th>lineage</th>"+
    "<th>Jaccard</th><th>overlap bp</th><th></th></tr>"+body+"</table>"+
    '<div class="note">Jaccard and overlap are against this locus\u2019 hg38 span. '+
    'ERV Navigator coordinates were published on hg19 and lifted here.</div></div>';
}
function hml2Panel(rows){
  if(!rows||!rows.length)return "";
  return rows.map(r=>'<div class="panel"><h2>HML-2 provirus detail — '+
    esc(r.subramanian_id||"")+'</h2><dl class="kv">'+
    kv("estimated age (MYA)",r.estimated_age_mya)+
    kv("oldest common ancestor",r.oldest_common_ancestor)+
    kv("ORFs",r.orfs||"none reported")+
    kv("polymorphic",r.polymorphic)+
    kv("hg19 (as published)",r.chrom_hg19?r.chrom_hg19+":"+fmt(r.start_hg19)+"-"+fmt(r.end_hg19):null)+
    kv("strand",r.strand)+
    kv("mapping",r.map_status+(r.start_offset_bp!=null?" · start offset "+r.start_offset_bp+" bp":""))+
    kv("source table",r.source_table==="T1"?"Table 1 (fairly intact)":"Table 2 (partial internal)")+
    (r.is_tandem_partner?kv("note","one of a tandem pair sharing this locus"):"")+
    '</dl><div class="note">Subramanian et al. 2011, Retrovirology 8:90 — '+
    'coordinates lifted hg19→hg38 by two-way liftover consensus.</div></div>').join("");
}

function coordTable(d,xg){
  const rows=(d.coord||[]).map(c=>{
    const u=UCSC[c.assembly]; const pos=c.chrom+":"+(c.start+1)+"-"+c.end;
    return "<tr><td>"+esc(c.assembly)+'</td><td class=mono>'+
      (u?'<a target="_blank" href="'+u+pos+'">'+esc(pos)+"</a>":esc(pos))+
      "</td><td>"+esc(c.strand)+"</td><td class=mono>"+bp(c.span)+"</td><td>"+esc(c.contig_type||"—")+
      "</td><td>"+(c.coord_unreliable?"unreliable":"ok")+"</td><td>"+esc(c.lift_method||"—")+"</td></tr>";}).join("");
  return "<table><tr><th>assembly</th><th>position (1-based display)</th><th>strand</th><th>span</th>"+
    "<th>contig</th><th>flag</th><th>lift</th></tr>"+rows+"</table>"+
    (xg.pct_identity!=null?'<div class="note">cross-genome: '+Number(xg.pct_identity).toFixed(2)+
      "% identity, "+esc(xg.xg_class||"")+", edit distance "+fmt(xg.edit_dist)+"</div>":"")+
    '<div class="note">Stored coordinates are 0-based half-open (UCSC/BED). Displayed as 1-based inclusive.</div>';
}
function aliasTable(al){
  const by={}; al.forEach(a=>{(by[a.alias_type]=by[a.alias_type]||[]).push(a);});
  return '<div class="acols">'+Object.keys(by).sort().map(t=>{
    const seen=new Set(),out=[];
    by[t].forEach(a=>{const k=a.alias+"|"+(a.assignment||"");if(seen.has(k))return;seen.add(k);
      out.push('<span class="mono">'+esc(a.alias)+"</span>"+
        (a.assignment?' <span class="badge">'+esc(a.assignment)+"</span>":"")+
        (a.is_current?"":' <span class="badge retired">retired</span>'));});
    return '<div class="ai"><span>'+esc(t)+"</span><span>"+out.join("<br>")+"</span></div>";}).join("")+"</div>";
}

// ---- locus graphic ----
async function drawLocus(d,co,asm){
  // called unawaited from render(); a throw here would otherwise be an invisible
  // rejected promise, leaving the "rendering…" placeholder on screen forever.
  try{ await drawLocus_(d,co,asm); }
  catch(e){ console.error("drawLocus",e);
    const g=$("gfx"); if(g) g.innerHTML='<div class="note" style="color:#a33">graphic failed: '+
      esc(e&&e.message||String(e))+'</div>'; }
}
// RepeatMasker lane colouring. Every class is drawn. ERV elements take the SAME two
// colours as the segment lane (LTR vs internal), so an element keeps its colour from
// lane to lane; other TE classes are pale context tints; low-information repeats
// (simple, low complexity, satellite, RNA, unknown) are grey, half height, unlabelled.
const REP_TE={LINE:"rep-line","LINE?":"rep-line",SINE:"rep-sine","SINE?":"rep-sine",
  DNA:"rep-dna","DNA?":"rep-dna",Retroposon:"rep-other",RC:"rep-other"};
const isInternalName=n=>/(-int|_I|-I|_int)$/i.test(n||"");
function repRole(r){
  if(r.rep_class==="LTR"||r.rep_class==="LTR?") return isInternalName(r.rep_name)?"int":"ltr";
  return REP_TE[r.rep_class]||"rep-low";
}
const STRONG=new Set(["ltr","int"]);

async function drawLocus_(d,co,asm){
  asm=asm||"hg38"; LASTDRAW=[d,co,asm];
  const isHg=asm==="hg38";
  // L is the lane-label gutter. At L=62 only ~9 monospace chars fit and 81% of
  // gene labels were clipped (e.g. "LOC124905662" rendered as "OC124905662").
  // L=110 holds ~18; longer labels are ellipsised by laneLabel() with the full
  // string in a <title>, so nothing is silently truncated.
  const PAD=1000, W=1080, L=110, R=14;
  const w0=Math.max(0,co.start-PAD), w1=co.end+PAD, span=w1-w0;
  const x=p=>L+(Math.min(Math.max(p,w0),w1)-w0)/span*(W-L-R);
  // gEVE ORFs and HERVarium domains are stored in hg38 coordinates ONLY, so they
  // are omitted on t2t rather than drawn at wrong positions. Gene models are NOT:
  // the bundle is assembly-keyed, hg38 from GENCODE and t2t from hs1 RefSeq, each
  // in its own assembly's coordinates. Index by assembly, never assume hg38.
  const GA=((typeof LOOKUP!=="undefined"&&LOOKUP)||{}).gene_assemblies||[];
  const tx=((d.gm||{})[asm]||[]).filter(t=>t.txEnd>w0&&t.txStart<w1);
  const lanes=[], USED=new Set();
  const reps=(d.repeats||[]).filter(r=>r.assembly===asm&&r.end>w0&&r.start<w1);
  // locus_segment is Telescope-derived and hg38-only. On t2t the LTR-class rows from
  // locus_repeat stand in, which is why the lane label names its source.
  const segs=isHg?(d.segments||[]).filter(s=>s.end>w0&&s.start<w1):[];
  lanes.push({label:isHg?"segments":"locus extent",h:16,draw:()=>segs.length?segs.map(s=>{
      const c=s.segment_class==="ltr"?"var(--ltr)":"var(--int)"; USED.add(s.segment_class==="ltr"?"ltr":"int");
      const a=x(s.start),b=x(s.end),txt=(s.repName||"")+(s.segment_class==="ltr"?" (LTR)":"");
      const tp=card(s.repName||"segment",[["segment",(s.segment_class||"")+" #"+s.seg_index],
        ["family",[s.repFamily,s.repClass].filter(Boolean).join(" / ")],
        ["position",co.chrom+":"+ivx(s.start,s.end)+" ("+esc(s.strand||"")+")"],["length",bpx(s.end-s.start)],
        ["RepeatMasker SW score",s.rmsk_sw_score]],c);
      return rect(a,0,b-a,13,c,tp)+
        lbl((a+b)/2,9.5,txt,"var(--on-strong)",8,"middle",fits(b-a,txt,8));}).join("")
    :rect(x(co.start),0,x(co.end)-x(co.start),13,"var(--rep-low)",
       card("locus extent",[["position",co.chrom+":"+ivx(co.start,co.end)],["length",bpx(co.end-co.start)],
         ["note",isHg?"no segments stored for this locus":"segments are hg38-only"]]))+
     lbl((x(co.start)+x(co.end))/2,9.5,
         isHg?"locus extent (no segments stored)":"locus extent (t2t; segments are hg38-only)",
         "var(--on-rep)",8,"middle",true)});
  // RepeatMasker lane: both assemblies. Labels are length-aware and separated, same
  // rule as the domain lane -- a 300 bp Alu is ~4 px wide at this scale.
  if(reps.length)lanes.push({label:"RepeatMasker",h:15,draw:()=>{
      let last=-1e9;
      return reps.slice().sort((p,q)=>p.start-q.start).map(r=>{
        const a=x(r.start),b=x(r.end),txt=r.rep_name||"", role=repRole(r), low=role==="rep-low"; USED.add(role);
        const ok=!low&&fits(b-a,txt,7)&&a-last>3; if(ok)last=b;
        const tp=card(txt||"repeat",[["class / family",[r.rep_class,r.rep_family].filter(Boolean).join(" / ")],
          ["position",co.chrom+":"+ivx(r.start,r.end)+" ("+(r.strand||"")+")"],["length",bpx(r.end-r.start)],
          ["divergence",r.pct_div==null?null:r.pct_div+"% from consensus"],
          ["shared by",r.n_loci>1?r.n_loci+" catalog loci":null]],"var(--"+role+")");
        return rect(a,low?3:0,Math.max(1.5,b-a),low?6:12,"var(--"+role+")",tp)+
               lbl((a+b)/2,8.7,txt,STRONG.has(role)?"var(--on-strong)":"var(--on-rep)",7,"middle",ok);
      }).join("");}});
  // Mappability lanes. Two resources on hg38 (pm151 Panmask "easy" at 151bp, and
  // Umap k100 single-read uniqueness); on t2t only Umap, because Panmask has no
  // T2T release -- the lane is omitted there rather than drawn empty, which would
  // read as "unmappable" when it means "no such resource".
  //
  // Blocks are shipped pre-clipped to locus +/-1000, the same window drawn here.
  // Both interval sets mark GOOD regions, so a GAP is the unmappable state; the
  // lane draws a faint full-width track under the blocks to make gaps legible as
  // absence rather than as background.
  const MAPRS=[["pm151","pm151 easy (151b)","var(--map1)"],
               ["umap100","Umap k100 unique","var(--map2)"]];
  const mp=d.mappability||{};
  MAPRS.forEach(([rs,lab,col])=>{
    const key=rs+"_"+asm, st=(mp.stats||{})[key], bl=(mp.blocks||{})[key];
    if(!st) return;                       // resource not present for this assembly
    if(st.no_data){
      lanes.push({label:lab,h:12,draw:()=>rect(x(w0),3,x(w1)-x(w0),7,"var(--maptrack)")+
        lbl((x(w0)+x(w1))/2,9,"no data \u2014 contig not covered by this resource",
            "var(--mut)",7,"middle",true)});
      return;
    }
    const segs2=(bl||[]).filter(b=>b[1]>w0&&b[0]<w1); USED.add("map1");
    const pct=st.frac==null?"":" \u2014 "+(100*st.frac).toFixed(0)+"% of locus";
    lanes.push({label:lab,h:12,draw:()=>{
      const tp=card(lab,[["mappable",st.frac==null?null:(100*st.frac).toFixed(0)+"% of locus"],
        ["longest unmappable run",st.longest_unmap==null?null:bpx(st.longest_unmap)],
        ["blocks overlapping locus",st.blocks],
        ["5\u2032 50 bp",st.t5==null?null:(st.t5?"mappable":"NOT mappable")],
        ["3\u2032 50 bp",st.t3==null?null:(st.t3?"mappable":"NOT mappable")],
        ["reading","filled = uniquely mappable; gaps are unmappable"]],col);
      let s='<g'+tipA(tp)+'>'+rect(x(w0),4,x(w1)-x(w0),5,"var(--maptrack)");   // gaps show through as this
      for(const b of segs2){
        const a=x(b[0]),e2=x(b[1]);
        s+=rect(a,3,Math.max(0.8,e2-a),7,col);
      }
      // element extent markers, so a reader can tell locus-internal gaps from
      // flanking ones without cross-referencing another lane
      s+=line(x(co.start),1,x(co.start),11,"var(--ink)",0.6)+
         line(x(co.end),1,x(co.end),11,"var(--ink)",0.6)+"</g>";
      return s;}});
  });
  // FANTOM5 per-locus CAGE peaks, drawn as ticks in the 5 px gap above the segment
  // lane (no lane of their own: fewer than 3% of loci carry one). Filled = same
  // strand as the locus, open = antisense.
  const f5m=((d.fantom5||{}).marks||{})[asm]||[];
  const f5in=f5m.filter(m=>m[0]>=w0&&m[0]<w1);
  if(f5in.length){
    USED.add("cage");
    const segLane=lanes[0], prev=segLane.draw;
    segLane.draw=()=>prev()+f5in.map(m=>{
      const X=x(m[0]), s=m[1]===1;
      const tp=card("FANTOM5 CAGE peak",[["position",co.chrom+":"+(m[0]+1).toLocaleString()],
        ["strand",s?"same as locus (sense)":"opposite to locus (antisense)"],
        ["max TPM",m[2]],["libraries \u2265 1 TPM",m[3]],
        ["peak set",m[4]?"FANTOM5 hg38-reprocessing (new peak)":"FANTOM5 hg19, lifted"]],"var(--cage)");
      return '<path d="M'+(X-3.5)+' -5 L'+(X+3.5)+' -5 L'+X+' 0 Z" fill="'+(s?"var(--cage)":"var(--pan)")+
        '" stroke="var(--cage)" stroke-width="1"'+tipA(tp)+'/>';}).join("");
  }
  // FANTOM5 CAGE lane REMOVED at v0.8. The layer's hit rate was too low to
  // justify the vertical space: it drew on a small minority of loci while
  // costing 22px of graphic height on every one. It was also never
  // independent evidence -- FANTOM CAT clusters derive from the same primary
  // CAGE data, so the two lanes agreeing was one observation, not two.
  // The builder retains _fantom5_payload and --f5-parquet, and locus_fantom5
  // is untouched in the catalog, so restoring this lane needs no recomputation.
  const orfs=isHg?(d.geve||[]).filter(o=>o.hg38_orf_end>w0&&o.hg38_orf_start<w1):[];
  if(orfs.length)USED.add("orf");
  if(orfs.length)lanes.push({label:"gEVE ORFs",h:15,draw:()=>orfs.map(o=>{
      const a=x(o.hg38_orf_start),b=x(o.hg38_orf_end),txt=o.orf_class||"ORF";
      const tp=card(o.geve_orf_id||txt,[["ORF class",o.orf_class],
        ["position",o.hg38_chrom+":"+ivx(o.hg38_orf_start,o.hg38_orf_end)+" ("+(o.orf_strand||"")+")"],
        ["length",bpx(o.hg38_orf_end-o.hg38_orf_start)+" \u00b7 "+Math.floor((o.hg38_orf_end-o.hg38_orf_start)/3)+" codons"],
        ...covOrfRows(d,co,o.hg38_orf_start,o.hg38_orf_end)],"var(--orf)");
      return arrow(a,b,12,"var(--orf)",o.orf_strand,tp)+
        lbl((a+b)/2,9,txt,"var(--on-strong)",8,"middle",fits(b-a,txt,8));}).join("")});
  const doms=isHg?(d.domains||[]).filter(o=>o.hg38_end>w0&&o.hg38_start<w1):[];
  if(doms.length)USED.add("dom");
  if(doms.length)lanes.push({label:"HERVarium domains",h:15,draw:()=>{
      // labels only where the box holds them AND no drawn label is within 3px
      let last=-1e9;
      return doms.slice().sort((p,q)=>p.hg38_start-q.hg38_start).map(o=>{
        const a=x(o.hg38_start),b=x(o.hg38_end),txt=o.gene||o.domain_desc||"";
        const ok=fits(b-a,txt,8)&&a-last>3; if(ok)last=b;
        const tp=card((o.gene||"domain")+(o.domain_desc?" \u2014 "+o.domain_desc:""),[["element",o.element],
          ["status",o.status],["score",o.domain_score],
          ["position","hg38:"+ivx(o.hg38_start,o.hg38_end)],["length",bpx(o.hg38_end-o.hg38_start)]],"var(--dom)");
        return rect(a,0,Math.max(2,b-a),12,"var(--dom)",tp)+lbl((a+b)/2,8.7,txt,"var(--on-strong)",8,"middle",ok);
      }).join("");}});
  // HERVarium LTRs with U3 / R / U5 (GRCh38 only). U3 is the transcript 5' part, so it sits on
  // the right on the minus strand; the segment coordinates are HERVarium's own.
  const hl=isHg?(d.hervarium||[]).filter(e=>e.kind==="LTR"&&e.end>w0&&e.start<w1):[];
  if(hl.length){ USED.add("u3"); USED.add("rr"); USED.add("u5");
    lanes.push({label:"HERVarium LTRs",h:15,draw:()=>hl.map(e=>{
      const segs=[["U3",e.u3_start,e.u3_end,"var(--u3)"],["R",e.r_start,e.r_end,"var(--rr)"],["U5",e.u5_start,e.u5_end,"var(--u5)"]]
        .filter(z=>z[1]!=null&&z[2]!=null);
      const L_=n=>n==null?null:bpx(n);
      const tp=card(e.id,[["role",e.role],["U3 / R / U5 call",e.st+(e.conf!=null?" (confidence "+Number(e.conf).toFixed(2)+")":"")],
        ["U3",e.u3_start!=null?L_(e.u3_end-e.u3_start):null],["R",e.r_start!=null?L_(e.r_end-e.r_start):null],["U5",e.u5_start!=null?L_(e.u5_end-e.u5_start):null],
        ["PBS",e.pbs_start!=null?"at "+(e.pbs_start+1).toLocaleString():null],["PPT",e.ppt_start!=null?"at "+(e.ppt_start+1).toLocaleString():null],
        ["subfamily",e.sub],["position",e.chrom+":"+ivx(e.start,e.end)+" ("+e.strand+")"],
        ["locus link",e.via==="linked_internal"?"via its linked internal element":null]],"var(--rr)");
      const a=x(Math.max(e.start,w0)),b=x(Math.min(e.end,w1));
      let g='<g'+tipA(tp)+'>'+rect(a,1,Math.max(1,b-a),10,"transparent")+
        '<rect x="'+a+'" y="1" width="'+Math.max(1,b-a)+'" height="10" fill="none" stroke="var(--rr)" stroke-width="0.8"'+
        (e.st==="LOW_CONF"?' stroke-dasharray="2,1.5"':"")+'/>';
      segs.forEach(z=>{const sa=x(Math.max(z[1],w0)),sb=x(Math.min(z[2],w1));
        if(sb>sa) g+='<rect x="'+sa+'" y="1" width="'+(sb-sa)+'" height="10" fill="'+z[3]+'"/>'+lbl((sa+sb)/2,9,z[0],z[0]==="R"?"var(--on-strong)":"#1b1b1f",7.5,"middle",fits(sb-sa,z[0],7.5));});
      [e.pbs_start,e.ppt_start].forEach((p_,i)=>{ if(p_!=null&&p_>=w0&&p_<=w1){const X=x(p_);
        g+='<path d="M'+(X-3)+' 13.5 L'+(X+3)+' 13.5 L'+X+' 10.5 Z" fill="#1b1b1f"><title>'+(i?"PPT":"PBS")+'</title></path>';}});
      return g+"</g>";}).join("")});
  }
  // HERVOminer ORFs >= 81 aa, one thin row per strand x frame. Solid = identical to a gEVE ORF.
  const hp=isHg?"hg38_":"t2t_";
  const hv=(d.hvorf||[]).filter(o=>o[hp+"start"]!=null&&o[hp+"chrom"]===co.chrom&&o[hp+"end"]>w0&&o[hp+"start"]<w1);
  if(hv.length){ USED.add("hvo");
    lanes.push({label:"HERVOminer ORFs \u2265 81 aa",h:41,draw:()=>hv.map(o=>{
      const s0=o[hp+"start"],s1=o[hp+"end"],st=o[hp+"strand"],fr=st==="+"?s0%3:s1%3,row=(st==="+"?0:3)+fr,y0=row*7;
      const a=x(Math.max(s0,w0)),b=x(Math.min(s1,w1));
      const tp=card(o.orf_id,[["length",o.aa_len+" aa"],["position",co.chrom+":"+ivx(s0,s1)+" ("+st+")"],
        ["placed",o[hp+"placement"]==="exact"?"exact sequence match":"most similar stretch (identity "+Number(o[hp+"identity"]).toFixed(3)+")"],
        ["gEVE",isHg?(o.geve_orf_id||"no identical gEVE ORF"):null],...(isHg?covOrfRows(d,co,s0,s1):[])],"var(--hvo)");
      const solid=isHg&&o.geve_orf_id;
      return '<rect x="'+a+'" y="'+y0+'" width="'+Math.max(1,b-a)+'" height="5" fill="'+(solid?"var(--hvo)":"var(--pan)")+
        '" stroke="var(--hvo)" stroke-width="0.8"'+(o[hp+"placement"]==="similar"?' stroke-dasharray="2,1.2"':"")+tipA(tp)+'/>';}).join("")});
  }
  // collapse isoforms to one lane per gene: union of exons, widest tx extent.
  // Drawing 6 lanes of the same gene wastes vertical space and hides whether any exon is in view.
  const byGene=new Map();
  tx.forEach(t=>{const k=t.name2||t.name;
    if(!byGene.has(k))byGene.set(k,{k,n:0,a:t.txStart,b:t.txEnd,strand:t.strand,ex:[]});
    const g=byGene.get(k); g.n++; g.a=Math.min(g.a,t.txStart); g.b=Math.max(g.b,t.txEnd);
    const es=t.exonStarts||[],ee=t.exonEnds||[];
    for(let i=0;i<es.length;i++) if(ee[i]>w0&&es[i]<w1) g.ex.push([es[i],ee[i]]);});
  [...byGene.values()].sort((p,q)=>p.a-q.a).forEach(g=>{
    g.ex.sort((p,q)=>p[0]-q[0]);
    const mg=[]; for(const e of g.ex){const l=mg[mg.length-1];
      if(l&&e[0]<=l[1])l[1]=Math.max(l[1],e[1]);else mg.push([e[0],e[1]]);}
    USED.add("gene");
    const gov=(d.genes||[]).find(h=>h.genome===asm&&(h.ref_gene_name===g.k||h.ref_gene_id===g.k))||{};
    const unnamed=/^ENSG\d+/.test(g.k);
    const GT={protein_coding:"coding",lncRNA:"lncRNA",processed_pseudogene:"pseudogene",
      unprocessed_pseudogene:"pseudogene",transcribed_unprocessed_pseudogene:"pseudogene"};
    const glab=unnamed?"unnamed "+(GT[gov.gene_type]||"gene"):g.k+(g.n>1?" ("+g.n+" tx)":"");
    lanes.push({label:glab,h:15,draw:()=>{
      const ov=(d.genes||[]).find(h=>h.genome===asm&&(h.ref_gene_name===g.k||h.ref_gene_id===g.k))||{};
      const tp=card(g.k,[["transcripts in window",g.n],["strand",g.strand],
        ["span",co.chrom+":"+ivx(g.a,g.b)],["exons in window",mg.length||"none (intronic)"],
        ["type",ov.gene_type],["overlap with locus",ov.overlap_bp==null?null:bpx(ov.overlap_bp)+
          (ov.exon_overlap_bp?" ("+bpx(ov.exon_overlap_bp)+" exonic)":"")],
        ["source",isHg?"GENCODE v50":"RefSeq (hs1)"]],"var(--gene)");
      let s='<g'+tipA(tp)+'>'+rect(x(g.a),0,Math.max(1,x(g.b)-x(g.a)),12,"transparent")+
        line(x(g.a),6,x(g.b),6,"var(--gene)",1);
      for(const [a,b] of mg) s+=rect(x(a),1.5,Math.max(1.5,x(b)-x(a)),9,"var(--gene)");
      s+="</g>"+lbl(x(Math.min(g.b,w1))+4,9.5,g.strand,"var(--mut)",9,"start",true);
      // y=9.5 put this ON the gene line at y=6, which struck through the text.
      // 13 clears the line and still sits inside the 15px lane.
      if(!mg.length) s+=lbl((x(Math.max(g.a,w0))+x(Math.min(g.b,w1)))/2,13,
          "intron only \u2014 no exon in window","var(--gene)",7.5,"middle",true);
      return s;}});});
  // Snaptron arc lane, hg38 only (srav3h has no t2t build). Unshifted so arcs sit
  // above the feature lanes they span, matching the TU view's layout.
  // Two possible sources, and only one is present per build (--debug-local decides).
  // `arcs` = full records with exact sample counts; `pjx` = the bit-packed
  // reference, which decodes to the same row shape. Debug is preferred when both
  // somehow appear (a stale mixed bundle) because its values are exact.
  //
  // NOTE the two do not select the same junctions: `arcs` takes a flat top-CAP by
  // support across the whole window, `pjx` takes top-15 with an end inside the
  // element plus the strongest spanning junction on each strand. The packed set is the
  // element-relevant one -- it keeps in-element junctions the flat cap drops -- so a
  // debug page and a public page legitimately show different arcs. The lane caption
  // states which rule produced what is on screen.
  let _arcsrc=null;
  if(isHg){
    if(d.arcs&&d.arcs.jx&&d.arcs.jx.length) _arcsrc=d.arcs;
    else if(d.pjx&&d.pjx.w&&d.pjx.w.length&&typeof pjxDecode==="function")
      _arcsrc=pjxDecode(d.pjx,co);
  }
  if(_arcsrc&&_arcsrc.jx.length&&typeof arcLane==="function"){
    const AH=46, inwin=_arcsrc.jx.filter(j=>j[1]>w0&&j[0]<w1);
    const nAnti=inwin.filter(j=>j.length>7?!!j[7]:(co.strand==="+"||co.strand==="-")&&j[3]!==co.strand).length;
    const HB=nAnti?30:0;
    USED.add("jx");
    lanes.unshift({label:"junctions",sub:nAnti?"sense \u2191  antisense \u2193":"",h:AH+HB,
      draw:()=>line(L,AH,W-R,AH,"var(--line)",0.8)+arcLane(_arcsrc,x,w0,w1,AH,co.strand,HB,(d.ts||{}).jb)+
      (typeof tssMarks==="function"?tssMarks((d.tx||{})[asm],co,x,co.strand):"")});
  }
  // CCLE 31-mer segments: each locus-unique segment, shaded by the share of the
  // 1,019 CCLE cell lines in which its k-mers are found (count >= 2)
  const _ccm=(LOOKUP||{}).ccle||{}, cc=d.cc||{}, NL=_ccm.n_lines_kmer||_ccm.n_lines||1019;
  if(cc.seg&&cc.seg.length&&cc.asm===asm){
    USED.add("ccle");
    const s0=co.start, span0=co.end-co.start;
    const side=r=>r<0?(co.strand==="-"?"3\u2032 flank":co.strand==="+"?"5\u2032 flank":"left flank"):
                   r>=span0?(co.strand==="-"?"5\u2032 flank":co.strand==="+"?"3\u2032 flank":"right flank"):"body";
    lanes.push({label:"CCLE 31-mer",h:12,draw:()=>rect(x(w0),5,x(w1)-x(w0),2,"var(--maptrack)")+
      cc.seg.filter(g=>s0+g[0]+g[1]>w0&&s0+g[0]<w1).map(g=>{
        const a=x(s0+g[0]), b=x(s0+g[0]+g[1]), f=g[2]/NL;
        const op=ccO(g[2],NL);
        const tp=card("locus-unique 31-mer segment",[["region",side(g[0])],
          ["position",co.chrom+":"+ivx(s0+g[0],s0+g[0]+g[1])],["length",bpx(g[1])],
          ["CCLE lines detecting",g[2].toLocaleString()+" of "+NL.toLocaleString()+" (k-mer count \u2265 2)"],
          ["max k-mer count",g[3]?Math.round(g[3]).toLocaleString():"0"]],"var(--ccle)");
        return '<g'+tipA(tp)+'>'+rect(a,1,Math.max(1.5,b-a),10,"#fff")+
          '<rect x="'+a+'" y="1" width="'+Math.max(1.5,b-a)+'" height="10" rx="1.5" fill="var(--ccle)" fill-opacity="'+
          op.toFixed(2)+'" stroke="var(--ccle)" stroke-width="0.6"/></g>';}).join("")});
  }
  // recount3 tissue breadth by region (5' background window | 5' flank | body | 3' flank |
  // 3' background window), one row per sample kind; darker = covered in more of that kind's
  // groups (tissue_map_layer.py, METHODS S22). Unique segments only, so stretches without a
  // unique segment have no measurement; per-base positions await the third cluster run.
  const _tm=((LOOKUP||{}).tissue||{}).map, _rg=(d.ts||{}).rg, _pm=(LOOKUP||{}).profile;
  if(isHg&&d.cov&&d.cv&&d.pf&&COV.mode==="cov"){ covLanes(d,co,lanes,USED,x,w0,w1,W,L,R); }
  else if(isHg&&_pm&&d.pf&&_tm){
    // 100-bp recount3 tissue profile (toolkit v3, tissue_profile_layer.py, D44): per bin, the
    // groups of each kind in which >= 10% of samples reach >= 0.1 per base per 1e9 aligned bases.
    // All reads, so bins in non-unique sequence (dashed) can carry multi-mapping reads.
    USED.add("cat1");
    const TG=(LOOKUP||{}).tissue.groups, pf=d.pf, n=pf[1], b00=co.start+pf[0];
    const unr=a_=>{const o=[];for(let j=0;j<a_.length;j+=2)for(let k=0;k<a_[j+1];k++)o.push(a_[j]);return o;};
    const TOPS=unr(pf[6]), MP=unr(pf[7]);
    const RL=["GTEx normal","TCGA normal","TCGA tumour","CCLE"];
    const KC=[0,1,2,3].map(k=>unr(pf[2+k]));
    _pm.kind_n.forEach((N,k)=>{
      if(!N) return;
      const col=k===3?"var(--ccle)":"var(--cat1)", C=KC[k];
      lanes.push({label:RL[k]+" (\u00a0"+N+")",h:10,draw:()=>{
        let g=rect(x(w0),4,x(w1)-x(w0),2,"var(--maptrack)"), i=0;
        while(i<n){
          let j=i; const nu=MP[i]<5;
          while(j+1<n&&C[j+1]===C[i]&&(MP[j+1]<5)===nu) j++;
          const a_=b00+100*i, e_=b00+100*(j+1);
          if(C[i]>0&&e_>w0&&a_<w1){
            let best=null; for(let t=i;t<=j;t++){const v=TOPS[t]; if(v!==-1&&(!best||v[1]>best[1])) best=v;}
            const mps=MP.slice(i,j+1), mpm=mps.reduce((p_,q_)=>p_+q_,0)/mps.length;
            const op=(0.15+0.85*C[i]/N)*(nu?0.5:1);
            const tp=card("recount3 coverage, "+RL[k]+" \u2014 100-bp bins",[["position",co.chrom+":"+ivx(a_,e_)+" ("+(j-i+1)+" bin"+(j>i?"s":"")+")"],
              ["groups covering it",C[i]+" of "+N+" "+RL[k]+" groups"],
              ["all kinds (first bin)",RL.map((l,q)=>l+" "+KC[q][i]+"/"+_pm.kind_n[q]).join(" \u00b7 ")],
              ["highest share (any kind)",best?TG[best[0]][1]+" "+(best[1]<10?"<1":Math.round(best[1]/10))+"%":"\u2014"],
              ["uniquely mappable (Umap k100)",Math.round(10*mpm)+"% of bases"+(nu?" \u2014 mostly not unique: may include multi-mapping reads":"")],
              ["rule","\u2265 "+Math.round(100*_pm.min_share)+"% of a group\u2019s samples at mean coverage \u2265 "+_pm.min_cov+" per base per 10\u2079 aligned bases (all reads)"]],col);
            const xa=x(Math.max(a_,w0)), xb=x(Math.min(e_,w1)), wd=Math.max(1,xb-xa);
            g+='<g'+tipA(tp)+'>'+rect(xa,0,wd,9,"#fff")+'<rect x="'+xa+'" y="0" width="'+wd+'" height="9" fill="'+col+'" fill-opacity="'+op.toFixed(2)+
              '"'+(nu?' stroke="'+col+'" stroke-width="0.6" stroke-dasharray="2,1.5"':'')+'/></g>';
          }
          i=j+1;
        }
        return g;}});
    });
  } else if(isHg&&_tm&&_rg){
    USED.add("cat1");
    const TG=(LOOKUP||{}).tissue.groups, s0=co.start;
    const PN={body:"body (unique segments)",up:"5\u2032 flank (unique segments)",dn:"3\u2032 flank (unique segments)",
      bg_up:"5\u2032 background window, 2\u201310 kb out",bg_dn:"3\u2032 background window, 2\u201310 kb out"};
    const RL=["GTEx normal","TCGA normal","TCGA tumour","CCLE"];
    _tm.kinds.forEach((kn,k)=>{
      const N=_tm.kind_n[k]; if(!N) return;
      const col=k===3?"var(--ccle)":"var(--cat1)";
      lanes.push({label:RL[k]+" ("+N+")",h:10,draw:()=>rect(x(w0),4,x(w1)-x(w0),2,"var(--maptrack)")+
        _tm.parts.filter(p=>_rg[p]).map(p=>{
          const r=_rg[p], R_=r[0], a=s0+R_[0], b=s0+R_[R_.length-1]; if(b<=w0||a>=w1) return "";
          const n=r[2+k], op=n?0.15+0.85*n/N:0;
          const top=r.slice(6).map(t=>TG[t[0]][1]+" "+(t[1]<10?"<1":Math.round(t[1]/10))+"%").join(", ");
          const tp=card(PN[p]+" \u2014 recount3 coverage",[["position",co.chrom+":"+ivx(a,b)],["unique bp measured",bpx(r[1])+" in "+(R_.length/2)+" stretch"+(R_.length>2?"es":"")],
            ["groups covering it",RL.map((l,i)=>l+" "+r[2+i]+"/"+_tm.kind_n[i]).join(" \u00b7 ")],
            ["highest shares",top||"none"],
            ["rule","\u2265 "+Math.round(100*_tm.min_share)+"% of a group\u2019s samples at mean coverage \u2265 "+_tm.min_cov+" per base per 10\u2079 aligned bases"]],col);
          let g='';
          for(let i=0;i<R_.length;i+=2){ const ra=s0+R_[i], rb=s0+R_[i+1]; if(rb<=w0||ra>=w1) continue;
            const xa=x(Math.max(ra,w0)), xb=x(Math.min(rb,w1)), wd=Math.max(1.5,xb-xa);
            g+=rect(xa,0,wd,9,"#fff")+'<rect x="'+xa+'" y="0" width="'+wd+'" height="9" rx="1" fill="'+col+'" fill-opacity="'+
              op.toFixed(2)+'" stroke="'+col+'" stroke-width="0.6"/>'; }
          return '<g'+tipA(tp)+'>'+g+'</g>';}).join("")});
    });
  }
  let y=0,body="";
  lanes.forEach(ln=>{
    body+='<g transform="translate(0,'+y+')">'+
      laneLabel(L-6,10,ln.label,"var(--mut)",9.5,L-6)+
      (ln.sub?lbl(L-6,23,ln.sub,"var(--mut)",8.5,"end",true):"")+ln.draw()+"</g>";
    y+=ln.h+5;});
  // locus extent guides + axis
  const guides=line(x(co.start),0,x(co.start),y,"var(--guide)",1,"2,2")+
               line(x(co.end),0,x(co.end),y,"var(--guide)",1,"2,2");
  let axis=line(L,y+4,W-R,y+4,"var(--axis)",1);
  const ticks=5;
  // edge ticks anchor inward: a centred label at the last tick overflows the viewBox
  // and is clipped by the browser (the right-hand coordinate showed as "4,043,77").
  for(let i=0;i<=ticks;i++){const p=w0+(w1-w0)*i/ticks;
    const an=i===0?"start":(i===ticks?"end":"middle");
    axis+=line(x(p),y+4,x(p),y+8,"var(--axis)",1)+
      lbl(x(p),y+19,Math.round(p).toLocaleString(),"var(--mut)",9,an,true);}
  axis+=lbl(L,y+33,asm+" "+co.chrom+"  ·  window "+(w1-w0).toLocaleString()+" bp  ·  locus "+
        (co.end-co.start).toLocaleString()+" bp","var(--mut)",9.5,"start",true);
  const LEG=[["ltr","ERV LTR"],["int","ERV internal"],["orf","gEVE ORF"],["dom","HERVarium domain"],
    ["u3","U3"],["rr","R"],["u5","U5 (HERVarium; dashed outline = low-confidence call)"],
    ["hvo","HERVOminer ORF \u2265 81 aa (solid = identical gEVE ORF; dashed = placed by similarity)"],
    ["gene","gene exon / intron"],["cage","FANTOM5 CAGE peak (open = antisense)"],["rep-line","LINE"],
    ["rep-sine","SINE"],["rep-dna","DNA"],["rep-other","other TE"],["rep-low","simple / low complexity"],
    ["map1","uniquely mappable"],["ccle","CCLE 31-mer segment (darker = detected in more cell lines)"],
    ["cat1","recount3 region coverage (darker = covered in more tissue / tumour groups)"]];
  $("gfx").dataset.used=[...USED].join(",");
  const legend='<div class="lg">'+LEG.filter(l=>USED.has(l[0])).map(l=>'<span><i style="background:var(--'+l[0]+')'+
    (l[0]==="rep-low"?";height:6px;vertical-align:1px":"")+'"></i>'+esc(l[1])+"</span>").join("")+
    (USED.has("ccle")?'<span class="jxkey">CCLE lines detecting <svg width="178" height="12">'+
      [[0,"0"],[1,"1"],[10,"10"],[100,"100"],[((LOOKUP||{}).ccle||{}).n_lines_kmer||1019,"all"]].map((p,i)=>'<rect x="'+(i*36+1)+'" y="1" width="14" height="10" rx="1.5" '+
        'fill="var(--ccle)" fill-opacity="'+ccO(p[0],((LOOKUP||{}).ccle||{}).n_lines_kmer||((LOOKUP||{}).ccle||{}).n_lines||1019).toFixed(2)+'" stroke="var(--ccle)" stroke-width="0.6"/>'+
        '<text x="'+(i*36+18)+'" y="10" font-size="9" fill="var(--mut)">'+p[1]+"</text>").join("")+"</svg></span>":"")+
    (USED.has("covheat")?'<span class="jxkey">mean coverage per base per 10\u2079 aligned bases '+covKey()+'</span>':"")+
    (USED.has("covhatch")?'<span class="jxkey"><svg width="22" height="12">'+COV_HATCH+'<rect x="1" y="1" width="20" height="10" fill="#de4968"/><rect x="1" y="1" width="20" height="10" fill="url(#covhatch)"/></svg> '+
      'under half the bin has unique 100-mers (Umap): reads there are shared with other copies, so dips and spikes follow the aligner, not expression</span>':"")+
    (USED.has("jx")?'<span class="jxkey">junction samples <svg width="150" height="12">'+
      [[10,"10"],[1000,"1k"],[100000,"100k"]].map((p,i)=>'<line x1="'+(i*50+2)+'" y1="6" x2="'+(i*50+22)+
        '" y2="6" stroke="var(--jx-within)" stroke-width="'+jxW(p[0]).toFixed(2)+'" opacity="'+jxO(p[0]).toFixed(2)+
        '" stroke-linecap="round"/><text x="'+(i*50+26)+'" y="10" font-size="9" fill="var(--mut)">'+p[1]+"</text>").join("")+
      '</svg> \u00b7 above = sense, below = antisense \u00b7 dashed = no mappable anchor</span>':"")+
    '<span style="margin-left:auto">hover for details \u00b7 click to pin and copy</span></div>';
  $("gfx").innerHTML=(isHg&&d.cov&&d.cv&&d.pf?covCtl(d):"")+'<svg viewBox="0 -6 '+W+" "+(y+48)+'" width="'+W+'">'+guides+body+axis+"</svg>"+legend+
    (tx.length?"":'<div class="note">no '+(isHg?"GENCODE":"RefSeq")+
      ' transcript in this window'+
      (GA.includes(asm)?"":" (no "+asm+" gene models in bundle)")+"</div>")+
    (d.segments.length?"":'<div class="note">This locus has no stored RepeatMasker segments — '+
      "only telescope-origin loci carry them. Bar shows the merged locus extent.</div>")+
    // Legend for the arc lane. The classes and the quantisation are only meaningful
    // for the packed source, so the caption states which source drew the lane and
    // how many junctions were kept out of how many were considered.
    (_arcsrc&&_arcsrc.jx.length?'<div class="note">junctions: '+
      '<span style="color:var(--jx-edge)">━</span> one end in element  '+
      '<span style="color:var(--jx-within)">━</span> both ends in element  '+
      '<span style="color:var(--jx-span)">━</span> element within intron  '+
      '· above the line = sense, below = antisense · width = samples, same scale on every locus '+
      '· dashed = neither 50 bp anchor mappable'+
      (_arcsrc.packed
        ? " · Snaptron srav3h, strand-aware, canonical only; top 15 with an end in the element, plus the strongest "+
          "spanning junction on each strand; support and depth log-quantised (\u00b14%)"
        : " · Snaptron srav3h, top "+_arcsrc.shown+" of "+
          _arcsrc.n_total.toLocaleString()+" by sample count, exact values")+
      "</div>":"");
}
const tipA=t=>t?' data-tip="'+esc(t)+'"':"";
const rect=(x,y,w,h,f,tip)=>'<rect x="'+x+'" y="'+y+'" width="'+Math.max(1,w)+'" height="'+h+'" fill="'+f+'" rx="1.5"'+tipA(tip)+'/>';
/* hover-card body. rows: [[key,value],...]; values are escaped here. */
function card(title,rows,sw){
  const cp=v=>'<button class="tcp" data-copy="'+esc(v)+'" title="copy">copy</button>';
  return (sw?'<i style="background:'+sw+'"></i>':"")+'<b class="tt">'+esc(title)+"</b>"+cp(title)+
    rows.filter(r=>r&&r[1]!=null&&r[1]!=="").map(r=>{
      const v=String(r[1]), pos=/^(position|span)$/.test(r[0]);
      // copy the bare interval, without the "(+)" strand suffix
      return '<br><span class="k">'+esc(r[0])+":</span> "+esc(v)+(pos?cp(v.replace(/\s*\([+-]?\)\s*$/,"").replace(/\u2013/g,"-")):"");
    }).join("");
}
const bpx=n=>Number(n).toLocaleString()+" bp";
const ivx=(a,b)=>(a+1).toLocaleString()+"\u2013"+Number(b).toLocaleString();
const line=(x1,y1,x2,y2,c,w,dash)=>'<line x1="'+x1+'" y1="'+y1+'" x2="'+x2+'" y2="'+y2+
  '" stroke="'+c+'" stroke-width="'+w+'"'+(dash?' stroke-dasharray="'+dash+'"':"")+"/>";
// Lane label, right-anchored in the L-px gutter. Ellipsises rather than letting the
// browser clip at the viewBox edge, and keeps the full name in a <title> tooltip.
const laneLabel=(x,y,t,c,s,gutter)=>{
  const per=0.60*s, max=Math.max(3,Math.floor((gutter-4)/per)), full=String(t||"");
  const cut=full.length>max?full.slice(0,max-1)+"\u2026":full;
  return '<text x="'+x+'" y="'+y+'" fill="'+c+'" font-size="'+s+
    '" text-anchor="end" font-family="ui-monospace,Menlo,monospace">'+
    (cut===full?"":'<title>'+esc(full)+'</title>')+esc(cut)+"</text>";
};
const lbl=(x,y,t,c,s,a,show)=>show?'<text x="'+x+'" y="'+y+'" fill="'+c+'" font-size="'+s+
  '" text-anchor="'+a+'" font-family="ui-monospace,Menlo,monospace">'+esc(t)+"</text>":"";
// monospace at font-size s is ~0.60*s per char; require the box to hold the string with padding
const fits=(w,t,s)=>w>=(String(t||"").length*0.60*s)+6;
function arrow(x1,x2,h,f,strand,tip){
  const w=x2-x1,t=Math.min(7,Math.max(2,w*0.25));
  if(w<5)return rect(x1,0,w,h,f,tip);
  return strand==="-"
    ?'<path d="M'+(x1+t)+' 0 H'+x2+' V'+h+' H'+(x1+t)+' L'+x1+' '+h/2+' Z" fill="'+f+'"'+tipA(tip)+'/>'
    :'<path d="M'+x1+' 0 H'+(x2-t)+' L'+x2+' '+h/2+' L'+(x2-t)+' '+h+' H'+x1+' Z" fill="'+f+'"'+tipA(tip)+'/>';
}

function hervariumPanel(d){
  const E=d.hervarium||[]; if(!E.length) return "";
  const ni=E.filter(e=>e.kind==="internal").length, nl=E.length-ni;
  const rows=E.slice().sort((p,q)=>p.start-q.start).map(e=>({"HERVarium name":e.id,kind:e.kind,subfamily:e.sub,
    role:e.kind==="LTR"?e.role:(e.ltr5_id||e.ltr3_id?"LTRs: "+(e.ltr5_id?"5\u2032":"")+(e.ltr5_id&&e.ltr3_id?" + ":"")+(e.ltr3_id?"3\u2032":""):"no linked LTR"),
    "U3/R/U5":e.kind==="LTR"?e.st:"", domains:e.kind==="internal"?e.nd:"",
    position:e.chrom+":"+ivx(e.start,e.end)+" ("+e.strand+")"}));
  return '<div class="panel"><h2>HERVarium elements \u2014 '+ni+' internal, '+nl+' LTR</h2>'+
    '<div class="note" style="font-style:normal">HERVarium v5 names (searchable). Internal elements and their linked LTRs; '+
    'LTRs carry HERVarium\u2019s U3/R/U5 call (OK or LOW_CONF). GRCh38 only.</div>'+
    tbl(rows,["HERVarium name","kind","subfamily","role","U3/R/U5","domains","position"])+"</div>";
}
function hvorfPanel(d){
  const O=d.hvorf||[]; if(!O.length) return "";
  const g=O.filter(o=>o.geve_orf_id).length, sim=O.filter(o=>o.hg38_placement==="similar").length;
  const rows=O.slice().sort((p,q)=>q.aa_len-p.aa_len).map(o=>({ORF:o.orf_id,aa:o.aa_len,
    hg38:o.hg38_start!=null?o.hg38_chrom+":"+ivx(o.hg38_start,o.hg38_end)+" ("+o.hg38_strand+")"+(o.hg38_placement==="similar"?" ~":""):"not placed",
    hs1:o.t2t_start!=null?o.t2t_chrom+":"+ivx(o.t2t_start,o.t2t_end)+" ("+o.t2t_strand+")"+(o.t2t_placement==="similar"?" ~":""):"not placed",
    gEVE:o.geve_orf_id||""}));
  return '<div class="panel"><h2>HERVOminer ORFs \u2265 81 aa \u2014 '+O.length+'</h2>'+
    '<div class="note" style="font-style:normal">'+g+' identical to a gEVE ORF. Positions are re-derived from each protein sequence; '+
    '\u201c~\u201d = placed on the most similar stop-to-stop stretch ('+sim+' here). Shorter ORFs are in the catalog only.</div>'+
    tbl(rows,["ORF","aa","hg38","hs1","gEVE"])+"</div>";
}

// Expression of the GENCODE v50 genes at this locus (+/- 1 kb), recount3 coverage over each
// gene's exons (all reads), per tissue group (gene_expression_layer.py, D44).
function geneExprPanel(d){
  const gx=d.gx, m=(LOOKUP||{}).gene_expr, T=((LOOKUP||{}).tissue||{});
  if(!gx||!gx.length||!m||!T.groups) return "";
  const TG=T.groups, KN=(T.map||{}).kind_n||[55,24,46,2], RL=["GTEx","TCGA normal","TCGA tumour","CCLE"];
  const rows=gx.map(g=>{
    const rel=g[4]||"";
    const tag=(/own/.test(rel)?"element\u2019s own gene model":(/contains/.test(rel)?"host ("+rel.replace(/ ?contains| ?own/g,"")+")":rel))||"\u2014";
    return {gene:g[1], type:g[2], strand:g[3], relation:tag,
      "groups expressing (median \u2265 1)":RL.map((l,k)=>l+" "+g[5+k]+"/"+KN[k]).join(" \u00b7 "),
      "highest medians":(g[9]||[]).map(t=>TG[t[0]][1]+" "+(t[1]/10).toFixed(1)).join(", ")||"\u2014"};});
  return '<div class="panel"><h2>Expression of genes at this locus \u2014 '+gx.length+'</h2>'+
    '<div class="note" style="font-style:normal">recount3 coverage over each GENCODE v50 gene\u2019s exons (all reads, so HERV bases inside the gene count too), '+
    'as mean coverage per exon base per 10\u2079 aligned bases; median per group, release sample exclusions applied. Genes overlapping the locus \u00b1 1 kb.</div>'+
    tbl(rows,["gene","type","strand","relation","groups expressing (median \u2265 1)","highest medians"])+"</div>";
}

// RNA Atlas Telescope TPM (RetroelementDB; Russ & Iordanskiy 2025; RNA Atlas GSE138734), D45.
// One row per organ system within a sample type; bars = share of samples >= 1 TPM in the polyA
// (solid) and total-RNA (outlined) library; diamond = recount3 coverage detection share in the
// matched GTEx tissue group(s) (tissue rows only), i.e. our alignment-free-of-Telescope readout.
function rnaAtlasPanel(d){
  const m=(LOOKUP||{}).rna_atlas, T=((LOOKUP||{}).tissue||{});
  if(!m) return "";
  const head='<div class="panel"><h2>RNA Atlas \u2014 Telescope TPM, 295 samples \u00d7 polyA and total RNA</h2>'+
    '<div class="note" style="margin:-2px 0 6px">RetroelementDB quantification (Russ &amp; Iordanskiy, <i>Mobile DNA</i> 2025) of the RNA Atlas (GSE138734)</div>';
  const isTel=(d.aliases||[]).some(a=>a.alias_type==="telescope_id");
  if(!isTel) return head+'<div class="note">not a Telescope locus \u2014 the RNA Atlas quantification covers Telescope loci only.</div></div>';
  const L=m.libs, val={}; (d.ra||[]).forEach(v=>{val[v[0]]=v[1]/100;});
  if(!d.ra||!d.ra.length) return head+'<div class="note">no library at \u2265 '+m.min_tpm_stored+' TPM.</div></div>';
  const TY=[["tissue","Tissues"],["cancer cell line","Cancer cell lines"],["cell line","Other cell lines"],["cell type","Primary cell types"]];
  const rows={};   // type|system -> {samples: {name: [polyA, total]}}
  L.forEach((l,i)=>{ const k=l[2]+"|"+l[3]; (rows[k]=rows[k]||{t:l[2],s:l[3],sm:{}});
    const o=(rows[k].sm[l[1]]=rows[k].sm[l[1]]||[null,null]); o[l[4]]=val[i]||0; });
  const gd=(d.ts||{}).d, TG=T.groups||[];
  const fmt1=v=>v==null?"\u2014":v>=100?Math.round(v).toLocaleString():v>=10?v.toFixed(0):v>=1?v.toFixed(1):v>0?v.toFixed(2):"0";
  const RH=15, LW=150, BW=150, W_=LW+BW+86;
  const col=(types)=>{
    let y=0, out="";
    types.forEach(([ty,tl])=>{
      const rs=Object.values(rows).filter(r=>r.t===ty).sort((a,b)=>a.s.localeCompare(b.s)); if(!rs.length) return;
      out+='<text x="0" y="'+(y+11)+'" font-size="10" font-weight="600" fill="var(--ink)">'+esc(tl)+'</text>'; y+=16;
      rs.forEach(r=>{
        const sm=Object.entries(r.sm), N=sm.length;
        const na=sm.filter(([,v])=>v[0]!=null).length, nt=sm.filter(([,v])=>v[1]!=null).length;
        const pa=sm.filter(([,v])=>v[0]>=m.detect_tpm).length, pt=sm.filter(([,v])=>v[1]>=m.detect_tpm).length;
        let gx=null, gl=[];
        if(ty==="tissue"&&gd){ const sh=[]; sm.forEach(([nm])=>{ const g=m.gtex[nm]; if(g){ const v=g.reduce((p_,q_)=>p_+gd[q_],0)/g.length/1000; sh.push(v); gl.push(nm.replace(" tissue","")+" \u2192 "+g.map(q_=>TG[q_][1]+" "+Math.round(gd[q_]/10)+"%").join(", ")); } });
          if(sh.length) gx=sh.reduce((p_,q_)=>p_+q_,0)/sh.length; }
        const top=sm.filter(([,v])=>(v[0]||0)>0||(v[1]||0)>0).sort((a,b)=>Math.max(b[1][0]||0,b[1][1]||0)-Math.max(a[1][0]||0,a[1][1]||0)).slice(0,8)
          .map(([nm,v])=>nm+": "+fmt1(v[0])+" / "+fmt1(v[1]));
        const tp=card(r.s+" \u2014 "+tl.toLowerCase()+" ("+N+" sample"+(N>1?"s":"")+")",[
          ["\u2265 "+m.detect_tpm+" TPM, polyA",pa+" of "+na],["\u2265 "+m.detect_tpm+" TPM, total RNA",pt+" of "+nt],
          ["highest (polyA / total TPM)",top.join("; ")||"none \u2265 "+m.min_tpm_stored],
          ...(gl.length?[["matched GTEx coverage (detected)",gl.join("; ")]]:[])],"var(--cat1)");
        const xa=v=>LW+BW*v;
        out+='<g'+tipA(tp)+'><rect x="0" y="'+y+'" width="'+W_+'" height="'+RH+'" fill="transparent"/>'+
          '<text x="'+(LW-6)+'" y="'+(y+10.5)+'" font-size="9.5" text-anchor="end" fill="var(--mut)">'+esc(r.s)+' ('+N+')</text>'+
          '<rect x="'+LW+'" y="'+(y+2)+'" width="'+BW+'" height="11" fill="var(--maptrack)"/>'+
          (na?'<rect x="'+LW+'" y="'+(y+2)+'" width="'+(BW*pa/na).toFixed(1)+'" height="5" fill="var(--cat1)"/>':'')+
          (nt?'<rect x="'+LW+'" y="'+(y+8)+'" width="'+Math.max(0,BW*pt/nt-0.6).toFixed(1)+'" height="4.4" fill="#fff" stroke="var(--cat1)" stroke-width="0.9"/>':'')+
          (gx!=null?'<path d="M'+xa(gx).toFixed(1)+' '+(y+1)+' l3.5 6.5 l-3.5 6.5 l-3.5 -6.5 z" fill="var(--ink)" fill-opacity="0.75"/>':'')+
          '<text x="'+(LW+BW+6)+'" y="'+(y+10.5)+'" font-size="9" fill="var(--mut)">'+pa+'/'+na+' \u00b7 '+pt+'/'+nt+'</text></g>';
        y+=RH+1; });
      y+=6; });
    return '<svg class="rachart" viewBox="0 0 '+W_+' '+y+'" width="'+W_+'" height="'+y+'" font-family="inherit">'+out+'</svg>';
  };
  const key='<div class="note" style="font-style:normal;margin:4px 0">'+
    '<svg width="22" height="10"><rect x="0" y="2" width="22" height="5" fill="var(--cat1)"/></svg> share of samples \u2265 '+m.detect_tpm+' TPM, polyA library \u00b7 '+
    '<svg width="22" height="10"><rect x="0.5" y="2" width="21" height="5" fill="#fff" stroke="var(--cat1)"/></svg> the same, total-RNA library \u00b7 '+
    '<svg width="10" height="12"><path d="M5 0 l4 6 l-4 6 l-4 -6 z" fill="var(--ink)" fill-opacity="0.75"/></svg> recount3 coverage detection (cpb &gt; 1) in the matched GTEx tissue(s) \u00b7 '+
    'numbers: samples \u2265 '+m.detect_tpm+' TPM, polyA \u00b7 total. Total RNA keeps unspliced and intronic RNA, so intronic loci light up there far more than in polyA or in GTEx (polyA).</div>';
  return head+key+'<div style="display:flex;gap:18px;flex-wrap:wrap;align-items:flex-start">'+col(TY.slice(0,3))+col(TY.slice(3))+"</div></div>";
}

// ---- coverage view (S28, D46) ---------------------------------------------------------
// d.cov: uint8 codes, groups x bins (page group order), bins of 100 bp from co.start+pf[0];
// value = 0.01 (2^(c/12) - 1) mean coverage per base per 1e9 aligned bases (all reads).
let COV={mode:"det",sel:{}}, LASTDRAW=null;
const COV_KIND={"GTEx normal":0,"TCGA adjacent normal":1,"TCGA tumour":2,"TCGA metastatic":2,"cell line (CCLE)":3};
const COV_RAMP=["#fcfdbf","#fed395","#fea772","#f7705c","#de4968","#b73779","#8c2981","#641a80","#3b0f70","#140e36"];
const COV_LINE=["#D55E00","#0072B2","#009E73","#CC79A7","#E69F00","#56B4E9","#1b1b1f","#999999"];
const covV=c=>c?0.01*(Math.pow(2,c/12)-1):0;
const covF=v=>v==null?"\u2014":v>=100?Math.round(v).toLocaleString():v>=10?v.toFixed(0):v>=1?v.toFixed(1):v>=0.01?v.toFixed(2):"<0.01";
function covCol(v){ if(v<0.01) return null; const t=Math.min(0.999,(Math.log10(v)+2)/(Math.log10(50)+2)); return COV_RAMP[Math.floor(t*COV_RAMP.length)]; }
function covKey(){ return '<svg width="236" height="12">'+[0.01,0.03,0.1,0.3,1,3,10,30].map((v,i)=>'<rect x="'+(i*29+1)+'" y="1" width="12" height="10" fill="'+covCol(v*1.05)+'"/>'+
  '<text x="'+(i*29+15)+'" y="10" font-size="8.5" fill="var(--mut)">'+(v<1?String(v).replace("0.","."):v)+'</text>').join("")+"</svg>"; }
function covUnr(a){ const o=[]; for(let j=0;j<a.length;j+=2) for(let k=0;k<a[j+1];k++) o.push(a[j]); return o; }
const covAt=(d,g,i)=>covV(d.cov[g*d.cv[1]+i]);
function covMean(d,co,g,a,b){
  const b00=co.start+d.pf[0], n=d.cv[1]; let s=0,w=0;
  for(let i=Math.max(0,Math.floor((a-b00)/100)); i<n&&b00+100*i<b; i++){
    const lo=Math.max(a,b00+100*i), hi=Math.min(b,b00+100*i+100); if(hi>lo){ s+=covAt(d,g,i)*(hi-lo); w+=hi-lo; } }
  return w?s/w:null;
}
function covExonRest(d,co,g){
  const b00=co.start+d.pf[0], ex=covUnr(d.cv[2]); let se=0,ne=0,sr=0,nr=0;
  ex.forEach((f,i)=>{ const a=b00+100*i; if(a+100<=co.start||a>=co.end) return; const v=covAt(d,g,i); if(f){se+=v;ne++;} else {sr+=v;nr++;} });
  return (ne&&nr)?[se/ne,sr/nr]:null;
}
// share of each 100-bp bin covered by Umap k100 unique blocks (the map's mappability lane);
// null when the locus has no Umap data. Bins under 0.5 are hatched in the coverage view (D47).
function covUniq(d,co){
  const mp=(d.mappability||{}), bl=(mp.blocks||{}).umap100_hg38, stt=((mp.stats||{}).umap100_hg38)||{};
  if(!bl||stt.no_data) return null;
  const b00=co.start+d.pf[0], n=d.cv[1], out=new Float32Array(n);
  for(const [a,b] of bl){ for(let i=Math.max(0,Math.floor((a-b00)/100)); i<n&&b00+100*i<b; i++){
    const lo=Math.max(a,b00+100*i), hi=Math.min(b,b00+100*i+100); if(hi>lo) out[i]+=(hi-lo)/100; } }
  return out;
}
const COV_HATCH='<defs><pattern id="covhatch" width="4" height="4" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">'+
  '<rect width="4" height="4" fill="#ffffff" fill-opacity="0.35"/><line x1="0" y1="0" x2="0" y2="4" stroke="#3a3a3a" stroke-width="1.1" stroke-opacity="0.55"/></pattern></defs>';
function covSel(d){ return COV.sel[d.uid]||(COV.sel[d.uid]=d.cv[3].slice()); }
function covRedraw(){ if(LASTDRAW) drawLocus(...LASTDRAW); }
function covMode(m){ COV.mode=m; covRedraw(); }
function covAdd(g){ g=+g; if(!LASTDRAW||isNaN(g)) return; const s=covSel(LASTDRAW[0]); if(!s.includes(g)&&s.length<8) s.push(g); covRedraw(); }
function covDel(g){ if(!LASTDRAW) return; const s=covSel(LASTDRAW[0]); const i=s.indexOf(+g); if(i>=0) s.splice(i,1); covRedraw(); }
if(typeof window!=="undefined"){ window.covMode=covMode; window.covAdd=covAdd; window.covDel=covDel; }
function covCtl(d){
  const TG=(LOOKUP.tissue||{}).groups||[], on=COV.mode==="cov", sel=covSel(d);
  const b=(m,t)=>'<button class="covbtn'+(COV.mode===m?' on':'')+'" onclick="covMode(\''+m+'\')">'+t+'</button>';
  let h='<div class="covctl"><span class="k">tissue rows:</span>'+b("det","detection")+b("cov","coverage");
  if(on){
    h+='<span class="k" style="margin-left:14px">tracks:</span>'+sel.map((g,i)=>'<span class="covchip"><i style="background:'+COV_LINE[i%8]+'"></i>'+
      esc(TG[g][1])+' <a href="javascript:void(0)" onclick="covDel('+g+')" title="remove">\u00d7</a></span>').join("")+
      (sel.length<8?'<select onchange="covAdd(this.value)"><option value="">+ add group\u2026</option>'+
        TG.map((g,i)=>sel.includes(i)?"":'<option value="'+i+'">'+esc(g[1])+' ('+g[3]+')</option>').join("")+'</select>':"");
  }
  return h+"</div>";
}
function covOrfRows(d,co,a,b){
  if(!d.cov||!d.cv||!d.pf) return [];
  const TG=(LOOKUP.tissue||{}).groups||[], sel=covSel(d).slice(0,4);
  const v=sel.map(g=>[g,covMean(d,co,g,a,b),covMean(d,co,g,co.start,co.end)]);
  return [["mean coverage over ORF (locus body)",v.map(t=>TG[t[0]][1]+" "+covF(t[1])+" ("+covF(t[2])+")").join(" \u00b7 ")]];
}
function covLanes(d,co,lanes,USED,x,w0,w1,W,L,R){
  const TG=(LOOKUP.tissue||{}).groups||[], n=d.cv[1], b00=co.start+d.pf[0], G=TG.length;
  const bm=TG.map((_,g)=>covMean(d,co,g,co.start,co.end)||0);
  const kinds=[[0,"GTEx"],[1,"TCGA normal"],[2,"TCGA tumour"],[3,"CCLE"]];
  const rows=[]; kinds.forEach(([k])=>{ TG.map((g,i)=>i).filter(i=>(COV_KIND[TG[i][2]]??2)===k).sort((p,q)=>bm[q]-bm[p]).forEach(i=>rows.push([i,k])); });
  const RH=1.8, GAP=4, KY={}; let yy=0, pk=-1;
  rows.forEach(r=>{ if(r[1]!==pk){ if(pk>=0) yy+=GAP; KY[r[1]]=[yy]; pk=r[1]; } r.push(yy); yy+=RH; KY[r[1]][1]=yy; });
  USED.add("covheat");
  const i0=Math.max(0,Math.floor((w0-b00)/100)), i1=Math.min(n,Math.ceil((w1-b00)/100));
  const UQ=covUniq(d,co), lowU=i=>UQ!=null&&UQ[i]<0.5;
  let nLow=0, nBody=0; for(let i=0;i<n;i++){ const a=b00+100*i; if(a+100>co.start&&a<co.end){ nBody++; if(lowU(i)) nLow++; } }
  if(nLow) USED.add("covhatch");
  const hatchRuns=(y0,h)=>{ let o="", i=i0; while(i<i1){ if(!lowU(i)){ i++; continue; } let j=i; while(j+1<i1&&lowU(j+1)) j++;
      const xa=x(Math.max(b00+100*i,w0)), xb=x(Math.min(b00+100*(j+1),w1)); o+='<rect x="'+xa.toFixed(1)+'" y="'+y0.toFixed(1)+'" width="'+Math.max(0.6,xb-xa).toFixed(1)+'" height="'+h.toFixed(1)+'" fill="url(#covhatch)" pointer-events="none"/>'; i=j+1; }
    return o; };
  lanes.push({label:"coverage \u00b7 "+G+" groups",h:Math.ceil(yy)+2,draw:()=>{
    let g="";
    kinds.forEach(([k,nm])=>{ if(KY[k]) g+=lbl(L-6,Math.max(KY[k][0]+7,(KY[k][0]+KY[k][1])/2+3),nm,"var(--mut)",8.5,"end",true); });
    rows.forEach(([gi,k,y0])=>{
      let runs="", i=i0;
      while(i<i1){ const c=covCol(covAt(d,gi,i)); let j=i; while(j+1<i1&&covCol(covAt(d,gi,j+1))===c) j++;
        if(c){ const xa=x(Math.max(b00+100*i,w0)), xb=x(Math.min(b00+100*(j+1),w1)); runs+='<rect x="'+xa.toFixed(1)+'" y="'+y0.toFixed(1)+'" width="'+Math.max(0.6,xb-xa).toFixed(1)+'" height="'+RH+'" fill="'+c+'"/>'; }
        i=j+1; }
      let pv=0,pi=i0; for(let q=i0;q<i1;q++){ const v=covAt(d,gi,q); if(v>pv){pv=v;pi=q;} }
      const er=covExonRest(d,co,gi);
      const tp=card(TG[gi][1]+" ("+TG[gi][3]+" samples)",[["mean over locus body",covF(bm[gi])],
        ["peak in window",covF(pv)+" at "+co.chrom+":"+ivx(b00+100*pi,b00+100*pi+100)],
        ["body exon bins : other body bins",er?covF(er[0])+" : "+covF(er[1]):"locus body has no GENCODE v50 exon"],
        ["hatched bins",UQ==null?"no Umap data":nLow+" of "+nBody+" body bins have < 50% unique 100-mers; coverage there depends on where the aligner placed shared reads"],
        ["unit","mean coverage per base per 10\u2079 aligned bases, mean over the group\u2019s samples (all reads)"]],covCol(Math.max(bm[gi],0.011))||"var(--maptrack)");
      g+='<g'+tipA(tp)+'><rect x="'+L+'" y="'+y0.toFixed(1)+'" width="'+(W-L-R)+'" height="'+RH+'" fill="var(--maptrack)"/>'+runs+'</g>';
    });
    return (nLow?COV_HATCH:"")+g+hatchRuns(0,yy);}});
  // tracks
  const sel=covSel(d), TH=104, top=6, bot=TH-4;
  let mx=0.1; sel.forEach(g=>{ for(let q=i0;q<i1;q++) mx=Math.max(mx,covAt(d,g,q)); });
  const lmax=Math.ceil(Math.log10(mx)), lmin=-2, yv=v=>bot-(Math.log10(Math.max(v,0.01))-lmin)/(lmax-lmin)*(bot-top);
  lanes.push({label:"",h:4,draw:()=>""});
  lanes.push({label:"coverage tracks",sub:"log scale",h:TH,draw:()=>{
    let g="";
    if(nLow){ let i=i0; while(i<i1){ if(!lowU(i)){ i++; continue; } let j=i; while(j+1<i1&&lowU(j+1)) j++;
      const xa=x(Math.max(b00+100*i,w0)), xb=x(Math.min(b00+100*(j+1),w1)); g+='<rect x="'+xa.toFixed(1)+'" y="0" width="'+Math.max(0.6,xb-xa).toFixed(1)+'" height="'+TH+'" fill="url(#covhatch)" opacity="0.6"/>'; i=j+1; } }
    for(let e=lmin;e<=lmax;e++){ const yy_=yv(Math.pow(10,e)); g+=line(L,yy_,W-R,yy_,"var(--axis)",0.4,"2,3")+lbl(W-R-2,yy_-2,e<0?String(Math.pow(10,e)).replace("0.","."):String(Math.pow(10,e)),"var(--mut)",8,"end",true); }
    sel.forEach((gi,s)=>{ let p=""; for(let q=i0;q<i1;q++){ const xa=x(Math.max(b00+100*q,w0)), xb=x(Math.min(b00+100*(q+1),w1)), yy_=yv(covAt(d,gi,q));
        p+=(q===i0?"M":"L")+xa.toFixed(1)+" "+yy_.toFixed(1)+"L"+xb.toFixed(1)+" "+yy_.toFixed(1); }
      g+='<path d="'+p+'" fill="none" stroke="'+COV_LINE[s%8]+'" stroke-width="1.3"/>'; });
    for(let q=i0;q<i1;q++){ const xa=x(Math.max(b00+100*q,w0)), xb=x(Math.min(b00+100*(q+1),w1));
      const tp=card("coverage at "+co.chrom+":"+ivx(b00+100*q,b00+100*q+100),sel.map(gi=>[TG[gi][1],covF(covAt(d,gi,q))]).concat(
        UQ==null?[]:[["unique 100-mers (Umap)",Math.round(100*UQ[q])+"% of the bin"+(UQ[q]<0.5?" \u2014 hatched: coverage depends on where shared reads were placed":"")]]),"var(--maptrack)");
      g+='<rect x="'+xa.toFixed(1)+'" y="0" width="'+Math.max(0.6,xb-xa).toFixed(1)+'" height="'+TH+'" fill="transparent"'+tipA(tp)+'/>'; }
    return g;}});
}

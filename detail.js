// Detail rendering + locus graphic. Depends on globals from index.html.
const UCSC={hs1:"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hs1&position=",
  hg38:"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=",
            t2t:"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hs1&position="};
let GENES=null;
async function geneModels(){
  if(GENES===null){
    try{const r=await fetch("data/gene_models.json.gz");
      GENES=r.ok?JSON.parse(await new Response(r.body.pipeThrough(new DecompressionStream("gzip"))).text()):{};}
    catch(e){GENES={};}
  }
  return GENES;
}
function coordOf(d,asm){return (d.coord||[]).find(c=>c.assembly===asm);}

// ---- v0.1 transcriptional units -------------------------------------------
// The TU layer is a SEPARATE assignment from the locus layer: units are built
// from element walks, not from loci, so one locus can span several units and a
// unit can carry several loci. Absent for loci excluded from v0.1.
function tuPanel(d){
  const t=d.tu||{};
  if(!t.units||!t.units.length){
    const why = t.in_v01===false ? "this locus is not represented in the v0.1 unit set"
                                 : "no v0.1 unit assignment in this bundle";
    return '<div class="panel"><h2>Transcriptional unit (v0.1)</h2>'+
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
  return '<div class="panel"><h2>Transcriptional unit (v0.1)</h2>'+
    '<div class="note">A separate assignment from the locus layer: units are built from '+
    'element walks, so a locus may span several units and a unit may carry several loci. '+
    'Cite the unit id only alongside the locus versioned_id.</div>'+warn+
    '<table class="tbl" style="margin-top:8px"><thead><tr><th>unit</th><th>group call</th>'+
    '<th>verdict</th><th>provenance</th><th>assemblies</th><th>internal bp</th>'+
    '<th>evidence</th><th>member groups</th></tr></thead><tbody>'+rows+'</tbody></table></div>';
}

function render(d,h){
  if(!d){$("view").innerHTML='<div class="empty">locus not found in bundle</div>';return;}
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
     '<div style="margin-top:9px">'+links+'</div></div>'+
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
   tuPanel(d)+
   // transcription evidence sits next to the TU panel: both describe what is
   // transcribed here, and the split question reads them together
   (typeof txPanel==="function"?txPanel(d):"")+
   (typeof jxPanel==="function"?jxPanel(d):"")+
   '<div class="panel"><h2>Locus map ±1 kb ('+(gasm||"—")+')'+
     (gco&&gasm!=="hg38"?' <span class="note" style="font-weight:400">— hg38 coordinate absent;'+
       ' gene models, gEVE ORFs and HERVarium domains are hg38-only and are omitted</span>':"")+
     '</h2>'+
     (gco?'<label class="note" style="display:block;margin:-4px 0 6px">'+
       '<input type="checkbox" id="repall"> show all repeat classes'+
       ' (simple repeats, low complexity, satellites)</label>':"")+
     '<div id="gfx">'+
     (gco?'<div class="note">rendering…</div>'
        :'<div class="note">no hg38 or t2t coordinate — graphic unavailable</div>')+'</div>'+
     '<div class="lg"><span><i style="background:var(--ltr)"></i>LTR segment</span>'+
     '<span><i style="background:var(--int)"></i>internal segment</span>'+
     '<span><i style="background:var(--orf)"></i>gEVE ORF</span>'+
     '<span><i style="background:var(--dom)"></i>HERVarium domain</span>'+
     '<span><i style="background:#b5651d"></i>LINE</span>'+
     '<span><i style="background:#2e8b8b"></i>SINE/Alu</span>'+
     '<span><i style="background:#3b6ea5"></i>LTR (RepeatMasker)</span>'+
     '<span><i style="background:#6b7f3a"></i>DNA</span>'+
     '<span><i style="background:var(--gene)"></i>gene exon (thick) / intron (thin)</span>'+

     '</div></div>'+
   '<div class="panel"><h2>Coordinates</h2>'+coordTable(d,xg)+'</div>'+
   '<div class="panel"><h2>Aliases — '+d.aliases.length+' rows</h2>'+aliasTable(d.aliases)+'</div>'+
   hml2Panel(d.hml2_detail)+
   missillacPanel(d.missillac)+
   '<div class="panel"><h2>Dfam best alignment</h2>'+
     (dfb.consensus_name?'<dl class="kv">'+kv("consensus",dfb.consensus_name)+kv("accession",dfb.dfam_accession)+
       kv("% identity",dfb.pct_identity==null?null:Number(dfb.pct_identity).toFixed(1))+
       kv("consensus coverage",dfb.cons_cov==null?null:(100*dfb.cons_cov).toFixed(1)+"%")+
       kv("SW score",dfb.sw_score)+kv("quality",dfb.aln_quality)+'</dl>'
      :'<div class="note">not screened against Dfam consensus</div>')+'</div>'+
   '<div class="panel"><h2>gEVE ORFs — '+d.geve.length+'</h2>'+tbl(d.geve,
      ["geve_orf_id","orf_class","hg38_chrom","hg38_orf_start","hg38_orf_end","orf_strand"])+'</div>'+
   '<div class="panel"><h2>HERVarium domains — '+d.domains.length+'</h2>'+tbl(d.domains,
      ["gene","domain_desc","element","status","domain_score","hg38_start","hg38_end","strand"])+'</div>'+
   '<div class="panel"><h2>Overlapping genes — '+d.genes.length+' rows</h2>'+tbl(d.genes,
      ["genome","source","ref_gene_name","ref_gene_id","gene_type","overlap_type","overlap_bp","exon_overlap_bp","n_transcripts"])+'</div>'+
   '<div class="panel"><h2>Segments — '+d.segments.length+'</h2>'+tbl(d.segments,
      ["seg_index","segment_class","repName","repFamily","repClass","chrom","start","end","strand","span","rmsk_sw_score"])+'</div>';
  if(gco) drawLocus(d,gco,gasm);
  const rc=$("repall"); if(rc) rc.onchange=()=>{ if(gco) drawLocus(d,gco,gasm); };
}
function kv(k,v){return "<dt>"+esc(k)+"</dt><dd>"+fmt(v)+"</dd>";}

// Snaptron splice-evidence summary. A strict ladder: the FIRST tier that
// applies is reported, so "donor or acceptor internal (sense)" implies no
// single sense junction had both ends inside, and an antisense tier implies no
// sense junction qualified at all.
//
// Two things this line deliberately does NOT do:
//  - it does not collapse "no evidence" into "not assessable". srav3h is
//    hg38-only and does not cover every alt contig, so 3,710 loci cannot be
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
  return "<table><tr><th>type</th><th>aliases</th></tr>"+Object.keys(by).sort().map(t=>{
    const seen=new Set(),out=[];
    by[t].forEach(a=>{const k=a.alias+"|"+(a.assignment||"");if(seen.has(k))return;seen.add(k);
      out.push('<span class="mono">'+esc(a.alias)+"</span>"+
        (a.assignment?' <span class="badge">'+esc(a.assignment)+"</span>":"")+
        (a.is_current?"":' <span class="badge retired">retired</span>'));});
    return "<tr><td>"+esc(t)+"</td><td>"+out.join("<br>")+"</td></tr>";}).join("")+"</table>";
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
// Repeat classes drawn by default. Everything else (Simple_repeat, Low_complexity,
// Satellite, tRNA/rRNA/snRNA/srpRNA, Unknown) is IN the shard and drawn only when the
// "all repeat classes" box is ticked -- see REP_LOW_INFO and the #repall handler.
const REP_INTERSPERSED=new Set(["LINE","SINE","LTR","LTR?","DNA","DNA?","Retroposon","RC","SINE?","LINE?"]);
const REP_COL={LINE:"#b5651d",SINE:"#2e8b8b","SINE?":"#2e8b8b",LTR:"#3b6ea5","LTR?":"#3b6ea5",
  DNA:"#6b7f3a","DNA?":"#6b7f3a",Retroposon:"#8f4c7a",RC:"#8f4c7a"};
const repCol=c=>REP_COL[c]||"#999";

async function drawLocus_(d,co,asm){
  asm=asm||"hg38";
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
  const GMB=await geneModels(), gm=(GMB&&GMB[asm])||{}, gkey=d.uid;
  const tx=(gm[gkey]||[]).filter(t=>t.txEnd>w0&&t.txStart<w1);
  const lanes=[];
  const showAll=(typeof document!=="undefined"&&$("repall")&&$("repall").checked);
  const reps=(d.repeats||[]).filter(r=>r.assembly===asm&&r.end>w0&&r.start<w1
                                       &&(showAll||REP_INTERSPERSED.has(r.rep_class)));
  // locus_segment is Telescope-derived and hg38-only. On t2t the LTR-class rows from
  // locus_repeat stand in, which is why the lane label names its source.
  const segs=isHg?(d.segments||[]).filter(s=>s.end>w0&&s.start<w1):[];
  lanes.push({label:isHg?"segments":"locus extent",h:16,draw:()=>segs.length?segs.map(s=>{
      const c=s.segment_class==="ltr"?"var(--ltr)":"var(--int)";
      const a=x(s.start),b=x(s.end),txt=(s.repName||"")+(s.segment_class==="ltr"?" (LTR)":"");
      return rect(a,0,b-a,13,c)+
        lbl((a+b)/2,9.5,txt,"#fff",8,"middle",fits(b-a,txt,8));}).join("")
    :rect(x(co.start),0,x(co.end)-x(co.start),13,"#c8c8d0")+
     lbl((x(co.start)+x(co.end))/2,9.5,
         isHg?"locus extent (no segments stored)":"locus extent (t2t; segments are hg38-only)",
         "#555",8,"middle",true)});
  // RepeatMasker lane: both assemblies. Labels are length-aware and separated, same
  // rule as the domain lane -- a 300 bp Alu is ~4 px wide at this scale.
  if(reps.length)lanes.push({label:"RepeatMasker",h:15,draw:()=>{
      let last=-1e9;
      return reps.slice().sort((p,q)=>p.start-q.start).map(r=>{
        const a=x(r.start),b=x(r.end),txt=r.rep_name||"";
        const ok=fits(b-a,txt,7)&&a-last>3; if(ok)last=b;
        return rect(a,0,Math.max(1.5,b-a),12,repCol(r.rep_class))+
               lbl((a+b)/2,8.7,txt,"#fff",7,"middle",ok);
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
  const MAPRS=[["pm151","pm151 easy (151b)","var(--map1,#2f7d4f)"],
               ["umap100","Umap unique (k100)","var(--map2,#4a6fa5)"]];
  const mp=d.mappability||{};
  MAPRS.forEach(([rs,lab,col])=>{
    const key=rs+"_"+asm, st=(mp.stats||{})[key], bl=(mp.blocks||{})[key];
    if(!st) return;                       // resource not present for this assembly
    if(st.no_data){
      lanes.push({label:lab,h:12,draw:()=>rect(x(w0),3,x(w1)-x(w0),7,"#f0eef2")+
        lbl((x(w0)+x(w1))/2,9,"no data \u2014 contig not covered by this resource",
            "#8a7f95",7,"middle",true)});
      return;
    }
    const segs2=(bl||[]).filter(b=>b[1]>w0&&b[0]<w1);
    const pct=st.frac==null?"":" \u2014 "+(100*st.frac).toFixed(0)+"% of locus";
    lanes.push({label:lab,h:12,draw:()=>{
      let s=rect(x(w0),4,x(w1)-x(w0),5,"#eceaef");   // gaps show through as this
      for(const b of segs2){
        const a=x(b[0]),e2=x(b[1]);
        s+=rect(a,3,Math.max(0.8,e2-a),7,col);
      }
      // element extent markers, so a reader can tell locus-internal gaps from
      // flanking ones without cross-referencing another lane
      s+=line(x(co.start),1,x(co.start),11,"#333",0.6)+
         line(x(co.end),1,x(co.end),11,"#333",0.6);
      s+='<title>'+esc(lab)+pct+
         (st.longest_unmap!=null?"; longest unmappable run "+st.longest_unmap.toLocaleString()+" bp":"")+
         (st.blocks!=null?"; "+st.blocks+" block(s) overlapping locus":"")+
         (st.t5!=null?"; 5\u2032 50bp "+(st.t5?"mappable":"NOT mappable"):"")+
         (st.t3!=null?"; 3\u2032 50bp "+(st.t3?"mappable":"NOT mappable"):"")+
         '</title>';
      return s;}});
  });
  // FANTOM5 CAGE lane REMOVED at v0.8. The layer's hit rate was too low to
  // justify the vertical space: it drew on a small minority of loci while
  // costing 22px of graphic height on every one. It was also never
  // independent evidence -- FANTOM CAT clusters derive from the same primary
  // CAGE data, so the two lanes agreeing was one observation, not two.
  // The builder retains _fantom5_payload and --f5-parquet, and locus_fantom5
  // is untouched in the catalog, so restoring this lane needs no recomputation.
  const orfs=isHg?(d.geve||[]).filter(o=>o.hg38_orf_end>w0&&o.hg38_orf_start<w1):[];
  if(orfs.length)lanes.push({label:"gEVE ORFs",h:15,draw:()=>orfs.map(o=>{
      const a=x(o.hg38_orf_start),b=x(o.hg38_orf_end),txt=o.orf_class||"ORF";
      return arrow(a,b,12,"var(--orf)",o.orf_strand)+
        lbl((a+b)/2,9,txt,"#fff",8,"middle",fits(b-a,txt,8));}).join("")});
  const doms=isHg?(d.domains||[]).filter(o=>o.hg38_end>w0&&o.hg38_start<w1):[];
  if(doms.length)lanes.push({label:"HERVarium domains",h:15,draw:()=>{
      // labels only where the box holds them AND no drawn label is within 3px
      let last=-1e9;
      return doms.slice().sort((p,q)=>p.hg38_start-q.hg38_start).map(o=>{
        const a=x(o.hg38_start),b=x(o.hg38_end),txt=o.gene||o.domain_desc||"";
        const ok=fits(b-a,txt,8)&&a-last>3; if(ok)last=b;
        return rect(a,0,Math.max(2,b-a),12,"var(--dom)")+lbl((a+b)/2,8.7,txt,"#fff",8,"middle",ok);
      }).join("");}});
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
    lanes.push({label:g.k+(g.n>1?" ("+g.n+" tx)":""),h:15,draw:()=>{
      let s=line(x(g.a),6,x(g.b),6,"var(--gene)",1);
      for(const [a,b] of mg) s+=rect(x(a),1.5,Math.max(1.5,x(b)-x(a)),9,"var(--gene)");
      s+=lbl(x(Math.min(g.b,w1))+4,9.5,g.strand,"#666",9,"start",true);
      // y=9.5 put this ON the gene line at y=6, which struck through the text.
      // 13 clears the line and still sits inside the 15px lane.
      if(!mg.length) s+=lbl((x(Math.max(g.a,w0))+x(Math.min(g.b,w1)))/2,13,
          "intron only \u2014 no exon in window","#7a5c8f",7.5,"middle",true);
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
  // element plus the single best spanning junction. The packed set is the
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
    const AH=46;
    lanes.unshift({label:"junctions",h:AH,draw:()=>arcLane(_arcsrc,x,w0,w1,AH)+
      (typeof tssMarks==="function"?tssMarks((d.tx||{})[asm],co,x,co.strand):"")});
  }
  let y=0,body="";
  lanes.forEach(ln=>{
    body+='<g transform="translate(0,'+y+')">'+
      laneLabel(L-6,10,ln.label,"#555",9.5,L-6)+ln.draw()+"</g>";
    y+=ln.h+5;});
  // locus extent guides + axis
  const guides=line(x(co.start),0,x(co.start),y,"#b9c6d4",1,"2,2")+
               line(x(co.end),0,x(co.end),y,"#b9c6d4",1,"2,2");
  let axis=line(L,y+4,W-R,y+4,"#999",1);
  const ticks=5;
  // edge ticks anchor inward: a centred label at the last tick overflows the viewBox
  // and is clipped by the browser (the right-hand coordinate showed as "4,043,77").
  for(let i=0;i<=ticks;i++){const p=w0+(w1-w0)*i/ticks;
    const an=i===0?"start":(i===ticks?"end":"middle");
    axis+=line(x(p),y+4,x(p),y+8,"#999",1)+
      lbl(x(p),y+19,Math.round(p).toLocaleString(),"#666",9,an,true);}
  axis+=lbl(L,y+33,asm+" "+co.chrom+"  ·  window "+(w1-w0).toLocaleString()+" bp  ·  locus "+
        (co.end-co.start).toLocaleString()+" bp","#666",9.5,"start",true);
  $("gfx").innerHTML='<svg viewBox="0 0 '+W+" "+(y+42)+'" width="'+W+'">'+guides+body+axis+"</svg>"+
    (tx.length?"":'<div class="note">no '+(isHg?"GENCODE":"RefSeq")+
      ' transcript in this window'+
      (Object.keys(gm).length?"":" (no "+asm+" gene models in bundle)")+"</div>")+
    (d.segments.length?"":'<div class="note">This locus has no stored RepeatMasker segments — '+
      "only telescope-origin loci carry them. Bar shows the merged locus extent.</div>")+
    // Legend for the arc lane. The classes and the quantisation are only meaningful
    // for the packed source, so the caption states which source drew the lane and
    // how many junctions were kept out of how many were considered.
    (_arcsrc&&_arcsrc.jx.length?'<div class="note">junctions: '+
      '<span style="color:'+JX_EDGE+'">━</span> one end in element  '+
      '<span style="color:'+JX_WITHIN+'">━</span> both ends in element  '+
      '<span style="color:'+JX_SPAN+'">━</span> element within intron  '+
      '· dashed = antisense to locus'+
      (_arcsrc.packed
        ? " · Snaptron srav3h, strand-aware, canonical only; top 15 in-element + top "+
          "spanning junction by sample count; support and depth log-quantised (\u00b14%)"
        : " · Snaptron srav3h, top "+_arcsrc.shown+" of "+
          _arcsrc.n_total.toLocaleString()+" by sample count, exact values")+
      "</div>":"");
}
const rect=(x,y,w,h,f)=>'<rect x="'+x+'" y="'+y+'" width="'+Math.max(1,w)+'" height="'+h+'" fill="'+f+'" rx="1.5"/>';
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
function arrow(x1,x2,h,f,strand){
  const w=x2-x1,t=Math.min(7,Math.max(2,w*0.25));
  if(w<5)return rect(x1,0,w,h,f);
  return strand==="-"
    ?'<path d="M'+(x1+t)+' 0 H'+x2+' V'+h+' H'+(x1+t)+' L'+x1+' '+h/2+' Z" fill="'+f+'"/>'
    :'<path d="M'+x1+' 0 H'+(x2-t)+' L'+x2+' '+h/2+' L'+(x2-t)+' '+h+' H'+x1+' Z" fill="'+f+'"/>';
}

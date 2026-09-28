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
   (d.genes.length?'<div class="panel"><h2>Overlapping genes — '+d.genes.length+' rows</h2>'+tbl(d.genes,
      ["genome","source","ref_gene_name","ref_gene_id","gene_type","overlap_type","overlap_bp","exon_overlap_bp","n_transcripts"])+'</div>':"")+
   (d.segments.length?'<div class="panel"><h2>Segments — '+d.segments.length+'</h2>'+tbl(d.segments,
      ["seg_index","segment_class","repName","repFamily","repClass","chrom","start","end","strand","span","rmsk_sw_score"])+'</div>':"")+
   absentPanel(d,dfb);
  if(gco) drawLocus(d,gco,gasm);
  document.querySelectorAll("button.fa").forEach(b=>b.onclick=()=>fetchFasta(d,b.dataset.fa,$("fastanote")));
  document.querySelectorAll("button.dl").forEach(b=>b.onclick=()=>downloadMap(d,gasm,b.dataset.dl));
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
function cclePanel(d){
  const cc=d.cc||{}, meta=(LOOKUP||{}).ccle; if(!meta||(!cc.tel&&!cc.bf)) return "";
  const NL=meta.n_lines, pct=n=>n==null?"\u2014":n.toLocaleString()+" ("+(100*n/NL).toFixed(n&&100*n/NL<1?1:0)+"%)";
  let tel="", telChart='<div></div>', kmChart='<div></div>';
  if(cc.tel){
    const t=cc.tel, TIS=meta.tissues;
    const chart=tisChart(TIS,t.tis.map(v=>v[2]),i=>[["at \u2265 1 TPM",null],["median TPM",t.tis[i][0]],["max TPM",t.tis[i][1]]],"at \u2265 1 TPM");
    tel='<div><h3>Telescope (TPM, hg38)</h3><dl class="kv">'+kv("feature",t.id)+kv("max TPM",t.mx)+kv("median TPM",t.md)+
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
      kv("body detected in",pct(bf.n_body))+
      kv(f5+" flank detected in",bf.n_5==null?"no unique flank segment":pct(bf.n_5))+
      kv(f3+" flank detected in",bf.n_3==null?"no unique flank segment":pct(bf.n_3))+
      kv("body > 2\u00d7 "+f5+" flank",bf.dom_5==null?"\u2014":pct(bf.dom_5))+
      kv("body > 2\u00d7 "+f3+" flank",bf.dom_3==null?"\u2014":pct(bf.dom_3))+
      kv("body > 2\u00d7 both flanks",bf.dom_both==null?"\u2014":pct(bf.dom_both))+
      (bf.hs!=null?kv("high-specificity set","yes \u2014 spec_ds "+bf.hs):"")+"</dl>"+
      '<div class="note">Detected = the region\u2019s weighted sum of unique 31-mer counts exceeds '+meta.det_wsum+
      " (counts below "+meta.min_abund+" ignored; weight \u221d segment length). Against Telescope this call is "+
      "99.8% specific (lines Telescope scores 0) and 98% sensitive (lines at \u2265 1 TPM). "+'Body vs flank compares length-normalised weighted sums (weight \u221d segment length), so a long body '+
      "is not favoured over a 1 kb flank by size."+(bf.oriented?"":" No locus strand: flanks are genomic left/right.")+
      " Segment positions are on the map (CCLE 31-mer lane).</div>"+
      "</div>";
    if(bf.tis) kmChart='<div><div class="note" style="margin-top:2px">31-mer: share of each tissue\u2019s cell lines with the body '+
      "detected (weighted sum > "+meta.det_wsum+", as above)</div>"+tisChart(meta.tissues,bf.tis,null,"body detected")+"</div>";
  } else km='<div><h3>31-mer (locus-unique k-mers)</h3><div class="note">no locus-unique 31-mer in the body: '+
    "k-mer evidence cannot be attributed to this locus (see mappability)</div></div>";
  return '<div class="panel"><h2>CCLE expression \u2014 '+NL.toLocaleString()+' cancer cell lines</h2><div class="ccgrid">'+
    tel+km+(cc.tel||(bf&&bf.tis)?telChart+kmChart:"")+'</div><div class="note" style="margin-top:8px">Two independent readouts of the same public CCLE RNA-seq runs: '+
    "Telescope reassigns multi-mapping reads by EM; the 31-mer readout asks whether k-mers unique to this locus occur in "+
    "each run. Where they disagree, the number of unique k-mers the locus has (mappability lanes, CCLE 31-mer lane) "+
    "is usually why.</div></div>";
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
  const LEG={ltr:"ERV LTR",int:"ERV internal",orf:"gEVE ORF",dom:"HERVarium domain",gene:"gene exon / intron",
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
function absentPanel(d,dfb){
  const miss=[];
  if(!dfb.consensus_name) miss.push(["Dfam alignment","not screened against Dfam consensus"]);
  if(!d.geve.length) miss.push(["gEVE ORFs",""]);
  if(!d.domains.length) miss.push(["HERVarium domains",""]);
  if(!d.genes.length) miss.push(["overlapping genes",""]);
  if(!d.segments.length) miss.push(["stored segments","only Telescope-origin loci carry them"]);
  if(!miss.length) return "";
  return '<div class="panel"><h2>Not present for this locus</h2><div class="absent">'+
    miss.map(m=>"<b>"+esc(m[0])+"</b>"+(m[1]?" ("+esc(m[1])+")":"")).join(" \u00b7 ")+"</div></div>";
}

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
        ["length",bpx(o.hg38_orf_end-o.hg38_orf_start)+" \u00b7 "+Math.floor((o.hg38_orf_end-o.hg38_orf_start)/3)+" codons"]],"var(--orf)");
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
      draw:()=>line(L,AH,W-R,AH,"var(--line)",0.8)+arcLane(_arcsrc,x,w0,w1,AH,co.strand,HB)+
      (typeof tssMarks==="function"?tssMarks((d.tx||{})[asm],co,x,co.strand):"")});
  }
  // CCLE 31-mer segments: each locus-unique segment, shaded by the share of the
  // 1,019 CCLE cell lines in which its k-mers are found (count >= 2)
  const cc=d.cc||{}, NL=((LOOKUP||{}).ccle||{}).n_lines||1019;
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
    ["gene","gene exon / intron"],["cage","FANTOM5 CAGE peak (open = antisense)"],["rep-line","LINE"],
    ["rep-sine","SINE"],["rep-dna","DNA"],["rep-other","other TE"],["rep-low","simple / low complexity"],
    ["map1","uniquely mappable"],["ccle","CCLE 31-mer segment (darker = detected in more cell lines)"]];
  $("gfx").dataset.used=[...USED].join(",");
  const legend='<div class="lg">'+LEG.filter(l=>USED.has(l[0])).map(l=>'<span><i style="background:var(--'+l[0]+')'+
    (l[0]==="rep-low"?";height:6px;vertical-align:1px":"")+'"></i>'+esc(l[1])+"</span>").join("")+
    (USED.has("ccle")?'<span class="jxkey">CCLE lines detecting <svg width="178" height="12">'+
      [[0,"0"],[1,"1"],[10,"10"],[100,"100"],[1019,"all"]].map((p,i)=>'<rect x="'+(i*36+1)+'" y="1" width="14" height="10" rx="1.5" '+
        'fill="var(--ccle)" fill-opacity="'+ccO(p[0],((LOOKUP||{}).ccle||{}).n_lines||1019).toFixed(2)+'" stroke="var(--ccle)" stroke-width="0.6"/>'+
        '<text x="'+(i*36+18)+'" y="10" font-size="9" fill="var(--mut)">'+p[1]+"</text>").join("")+"</svg></span>":"")+
    (USED.has("jx")?'<span class="jxkey">junction samples <svg width="150" height="12">'+
      [[10,"10"],[1000,"1k"],[100000,"100k"]].map((p,i)=>'<line x1="'+(i*50+2)+'" y1="6" x2="'+(i*50+22)+
        '" y2="6" stroke="var(--jx-within)" stroke-width="'+jxW(p[0]).toFixed(2)+'" opacity="'+jxO(p[0]).toFixed(2)+
        '" stroke-linecap="round"/><text x="'+(i*50+26)+'" y="10" font-size="9" fill="var(--mut)">'+p[1]+"</text>").join("")+
      '</svg> \u00b7 above = sense, below = antisense \u00b7 dashed = no mappable anchor</span>':"")+
    '<span style="margin-left:auto">hover for details \u00b7 click to pin and copy</span></div>';
  $("gfx").innerHTML='<svg viewBox="0 -6 '+W+" "+(y+48)+'" width="'+W+'">'+guides+body+axis+"</svg>"+legend+
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
